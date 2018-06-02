/**
 * proof of concept for using iopipe to parse fasta data
 *
 * Format:
 * >Entry1_ID header field1|header field2|...
 * CAGATATCTTTGATGTCCTGATTGGAAGGACCGTTGGCCCCCCACCCTTAGGCAG
 * TGTATACTCTTCCATAAACGAGCTATTAGTTATGAGGTCCGTAGATTGAAAAGGG
 * TGACGGAATTCGGCCGAACGGGAAAGACGGACATCTAGGTATCCTGAGCACGGTT
 * GCGCGTCCGTATCAAGCTCCTCTTTATAGGCCCCG
 * >Entry2_ID header field1|header field4|...
 * GTTACTGTTGGTCGTAGAGCCCAGAACGGGTTGGGCAGATGTACGACAATATCGCT
 * TAGTCACCCTTGGGCCACGGTCCGCTACCTTACAGGAATTGAGA
 *
 * >Entry3_ID header field1|header field2|...
 * GGCAGTACGATCGCACGCCCCACGTGAACGATTGGTAAACCCTGTGGCCTGTGAGC
 * GACAAAAGCTTTAATGGGAAATACGCGCCCATAACTTGGTGCGA
 *
 * Some characteristics:
 *
 * - Entry_ID is >[[:alphanumeric:]]. Where '>' marks the entry start.
 * - Headers may contain annotation information separated by some delimiter (i.e. | in this case).
 * - Entry ID and header is a single line, which does not contain newline characters.
 * - Sequence under the header line is [ATCGN\n]* (Perl regex).
 * - A fasta file can be plain-text or gzip compressed.
 */

/**
 * TODO
 * Implement
 * filter (i.e. by pattern on header),
 * writer
 * reverse complement
 *
 */


module fastaq.fasta;

import iopipe.traits;
import iopipe.textpipe;
private import std.traits;
private import std.range.primitives;
private import std.algorithm : find, splitter, filter;
private import std.conv: to;
private import std.string : stripLeft, stripRight, strip;

struct BufRef
{
    // position within the buffer of the starting reference
    size_t pos;
    // length of the reference
    size_t length;
    auto value(B)(B buf)
    {
        assert(pos <= buf.length);
        assert(pos + length <= buf.length);
        return buf[pos .. pos + length];
    }

    void release(size_t elements)
    {
        pos -= elements;
    }
}

struct FastaToken
{
    size_t endPos;
    BufRef entryid;
    BufRef[] fields;
    BufRef sequence;

    void release(size_t elements)
    {
        endPos -= elements;
        entryid.release(elements);
        sequence.release(elements);
        foreach(ref f; fields) f.release(elements);
    }

    auto value(B)(B buf)
    {
        FastaConcreteToken!B result;
        result.entryid = entryid.value(buf);
        result.fields = new B[fields.length];
        foreach(i, ref f; fields)
            result.fields[i] = f.value(buf);
        result.sequence = sequence.value(buf);
        return result;
    }
}

struct FastaConcreteToken(R)
{
    R entryid;
    R[] fields;
    R sequence;
}

auto tokenParser(Chain, char header = '>', char fieldsep = '|')(Chain c) if (isIopipe!Chain && isSomeChar!(ElementEncodingType!(WindowType!Chain)))
{
    auto lines = c.byLine;
    alias ChainType = typeof(lines);
    static struct Result
    {
        ChainType chain;
        size_t pos;
        alias chain this;

        FastaToken nextToken()
        {
            if(pos == chain.window.length)
                // reaches the end of the stream
                return FastaToken.init;
            // pos must start with a start identifier
            assert(chain.window[pos] == header);
            // the header is the current line
            FastaToken result;
            auto fields = chain.window[pos .. $].stripRight.splitter(fieldsep);
            if(!fields.empty)
            {
                auto firstElemSize = fields.front.length;
                auto firstField = fields.front.find(' ');
                result.entryid = BufRef(pos + 1, firstElemSize - firstField.length - 1);
                if(firstField.length > 0)
                {
                    firstField = firstField.stripLeft;
                    result.fields ~= BufRef(pos + (firstElemSize - firstField.length), firstField.length);
                }
                pos += firstElemSize;

                fields.popFront;
                pos += 1; // skip newline or |
                foreach(f; fields)
                {
                    if(!f.empty)
                        result.fields ~= BufRef(pos, f.length);
                    pos += f.length + 1;
                }
            }

            // parse all the sequence
            auto seqStart = pos;
            while(chain.extend(0) != 0)
            {
                if(chain.window[pos] == header)
                    break;
                pos = chain.window.length;
            }

            auto seqData = chain.window[seqStart .. pos].stripLeft;
            seqStart = pos - seqData.length;
            seqData = seqData.stripRight;
            result.sequence = BufRef(seqStart, seqData.length);
            result.endPos = pos;
            return result;
        }

        void release(size_t elements)
        {
            pos -= elements;
            chain.release(elements);
        }
    }

    // prime the lines item
    lines.extend(0);
    while (lines.window.strip.empty)
      {
        lines.release(lines.window.length);
        lines.extend(0);
      }
    return Result(lines);
}

unittest
{
    immutable auto input = "\n" ~ ">EntryId1 field1|field2|field3\n" ~
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" ~
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" ~
        "ACGTACGTACGTACGTACGTACG \n" ~
        "\n" ~
        ">EntryId2 field4|field5|length > 3\n" ~
        " ACGT \n" ~
        " ACG \n";

    auto tokenizer = input.tokenParser;
    auto item1 = tokenizer.nextToken;
    assert(item1.entryid.value(tokenizer.window) == "EntryId1");
    assert(item1.fields.length == 3);
    assert(item1.fields[0].value(tokenizer.window) == "field1");
    assert(item1.fields[1].value(tokenizer.window) == "field2");
    assert(item1.fields[2].value(tokenizer.window) == "field3");
    auto seq = item1.sequence.value(tokenizer.window);
    assert(seq[0] == 'A');
    assert(seq[$-1] == 'G');
    import std.range: cycle;
    import std.ascii: isWhite;
    import std.algorithm: filter, startsWith;
    assert(cycle("ACGT").startsWith(seq.filter!(a => !a.isWhite)));

    auto item2 = tokenizer.nextToken;

    assert(item2.entryid.value(tokenizer.window) == "EntryId2");
    assert(item2.fields.length == 3);
    assert(item2.fields[0].value(tokenizer.window) == "field4");
    auto field5 = item2.fields[1].value(tokenizer.window);
    assert(field5 == "field5", "got: " ~  field5);
    auto fieldspecial = item2.fields[2].value(tokenizer.window);
    assert(fieldspecial == "length > 3", "Expect 'length > 3' got: " ~ fieldspecial);
    seq = item2.sequence.value(tokenizer.window);
    assert(seq.filter!(a => !a.isWhite).to!string == "ACGTACG", "Expected: ACGTACG, got: " ~ seq);

    auto item3 = tokenizer.nextToken;
    assert(item3.entryid.length == 0);

    tokenizer.release(item1.endPos);
    item2.release(item1.endPos);

    auto concrete = item2.value(tokenizer.window);

    assert(concrete.entryid == "EntryId2");
    assert(concrete.fields.length == 3);
    assert(concrete.fields[0] == "field4");
    assert(concrete.fields[1] == "field5");
    seq = concrete.sequence;
    assert(seq.filter!(a => !a.isWhite).to!string == "ACGTACG", "Expected: ACGTACG, got: " ~ seq);
}
