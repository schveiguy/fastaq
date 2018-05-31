/**
 * High performance FASTQ library using iopipe
 *
 * Format:
 * Excerpt from https://en.wikipedia.org/wiki/FASTQ_format

 * FASTQ file normally uses four lines per sequence.
 * Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
 * Line 2 is the raw sequence letters.
 * Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
 * Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

 * A FASTQ file containing a single sequence might look like this:
 * @SEQ_ID
 * GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
 * +
 * !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

 * The character '!' represents the lowest quality while '~' is the highest.
 * Here are the quality value characters in left-to-right increasing order of
 * quality (ASCII):
 * !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

 *
 * Some characteristics:
 * - SEQ_ID is @[[:alphanumeric:]]. Where '@' marks the entry start. More description about Illumina header line below
 * - Entry ID and header is a single line, which does not contain newline characters.
 * - Sequence under the header line is [ATCGN]* (Perl regex).
 * - A fastq file can be plain-text or gzip compressed.

 * Illumina Header Line:
 * @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

 * EAS139	the unique instrument name
 * 136	the run id
 * FC706VJ	the flowcell id
 * 2	flowcell lane
 * 2104 	tile number within the flowcell lane
 * 15343 	'x'-coordinate of the cluster within the tile
 * 197393 	'y'-coordinate of the cluster within the tile
 * 1 	the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
 * Y 	Y if the read is filtered, N otherwise
 * 18 	0 when none of the control bits are on, otherwise it is an even number
 * ATCACG 	index sequence. It is sample number if it is a number
 *
 */

module fastaq.fastq;

import iopipe.traits;
import iopipe.textpipe;
private import std.traits;
private import std.range.primitives;
private import std.algorithm : find, splitter, filter;
private import std.conv: to;
private import std.ascii: isWhite;
private import std.string : stripLeft, stripRight, strip;
private import std.experimental.logger;

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

struct FastqToken
{
  size_t endPos;
  BufRef entryid;
  BufRef[] hfields;
  BufRef sequence;
  BufRef[] extras;
  BufRef optional;
  BufRef quality;

  void release(size_t elements)
  {
    endPos -= elements;
    entryid.release(elements);
    sequence.release(elements);
    optional.release(elements);
    quality.release(elements);
    foreach(ref f; hfields) f.release(elements);
    foreach(ref f; extras) f.release(elements);
  }

  auto value(B)(B buf)
  {
    FastqConcreteToken!B result;
    result.entryid = entryid.value(buf);
    result.hfields = new B[hfields.length];
    result.extras = new B[extras.length];
    foreach(i, ref f; hfields)
      result.hfields[i] = f.value(buf);
    foreach(j, ref e; extras)
      result.extras[j] = e.value(buf);
    result.sequence = sequence.value(buf);
    result.optional = optional.value(buf);
    result.quality = quality.value(buf);
    return result;
  }
}

struct FastqConcreteToken(R)
{
  R entryid;
  R[] hfields;
  R[] extras;
  R sequence;
  R optional;
  R quality;

}

auto tokenParser(Chain,
                 char fieldsep = ' ',
                 char subfs = ':')(Chain c) if (isIopipe!Chain && isSomeChar!(ElementEncodingType!(WindowType!Chain)))
  {
    auto lines = c.byLine;
    alias ChainType = typeof(lines);
    size_t lineCount = 0;
    static struct Result
    {
      ChainType chain;
      size_t pos;
      alias chain this;

      FastqToken nextToken()
      {
        // pos must start with a start identifier
        if(pos == chain.window.length)
          return FastqToken.init;
        // Empty line
        while (chain.window[0 .. $].strip.length == 0)
          {
            pos += chain.window.length;
            chain.extend(0);
          }
        logf("Pos: %s, window length: %s, window: %s", pos, chain.window.length, chain.window );
        assert(chain.window[pos] == '@', "Got this: " ~ chain.window[pos]);
        // the header is the current line
        pos += 1; // jump over the '@'
        FastqToken result;
        // parse header line
        auto lazyArr = chain.window[pos .. $].stripRight.splitter(fieldsep);
        pragma(msg, typeof(lazyArr).stringof);
        if(!lazyArr.empty)
          {
            auto firstElem = lazyArr.front;
            result.entryid = BufRef(pos, firstElem.length);
            auto i = 0;
            foreach(elem; lazyArr)
              {
                auto subFields = elem.splitter(subfs);
                foreach(f; subFields)
                  {
                    if (f.empty) continue;
                    if(i == 0)
                      {
                        result.hfields ~= BufRef(pos, f.length);
                      }
                    else
                      {
                        result.extras ~= BufRef(pos, f.length);
                      }
                    pos += f.length + 1;
                  }
                i += 1;
                // pos += 1; // skip fieldsep or newline at the end of line
              }
          }

        // parse all the sequence
        auto seqStart = pos;
        string seqData;
        if (chain.extend(0) != 0)
          {
            pos = chain.window.length;
            seqData = chain.window[seqStart .. pos].stripLeft;
            seqStart = pos - seqData.length;
            seqData = seqData.stripRight;
            result.sequence = BufRef(seqStart, seqData.length);
            // logf("Seqstart %s, pos: %s, seqData: %s", seqStart, pos, seqData);

          }
        // parse optional info line
        auto optLineStart = pos;
        string optLineData;
        if (chain.extend(0) != 0)
          {
            pos = chain.window.length;
            optLineData = chain.window[optLineStart .. pos].stripLeft;
            optLineStart = pos - optLineData.length;
            optLineData = optLineData.stripRight;
            result.optional = BufRef(optLineStart, optLineData.length);

            // logf("Opt line %s, pos: %s, content: %s", seqStart, pos, optLineData);
          }
        // parse quality
        auto qualLineStart = pos;
        string qualLineData;
        if (chain.extend(0) != 0)
          {
            pos = chain.window.length;
            qualLineData = chain.window[qualLineStart .. pos].stripLeft;
            qualLineStart = pos - qualLineData.length;
            qualLineData = qualLineData.stripRight;
            result.quality = BufRef(qualLineStart, qualLineData.length);
            // logf("Qual line %s, pos: %s, content: %s", qualLineStart, pos, qualLineData);

            result.endPos = pos;
            logf("Finish Rec, Pos: %s, window length: %s, window: %s", pos, chain.window.length, chain.window );
          }
        chain.extend(0);
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
    return Result(lines);
  }

version(testclass)
unittest{
  // Test make fastq entry and parse indexes
  bool test(){
    Read r;
    r = new Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                 "CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAG",
                 "+",
                 "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEE");
    string idx = r.lastIndex();
    writeln("first index: ", r.firstIndex);
    return idx == "GGTCCCGA" && r.firstIndex == "TATAGCCT";

  }

  assert(true == test());

}

/**
 * take a fastq entry and print out its content.
 *
 */
void printFastqEntry(T)(T fqentry){
  import std.stdio: writeln, writefln;
  writefln("Entry: %s", fqentry.stringof);
  writeln(fqentry.entryid);
  writeln(fqentry.hfields);
  writeln(fqentry.optional);
  writeln(fqentry.quality);
}
unittest
{
  auto input = r"
@M04916:6:000000000-BRH9H:1:1101:14943:1392 1:N:0:2
TTCTAATTCATCTTTGGAACAAGAACAGACCAAATGAGAAAAAATATATTTGAAGTTGTTTATTAAAAG
+
@@@B--CEE,CCEFF,,,,;,,,,,;,,,;C,,,C,,,,,,,,+;,<,CCE,,,,<C,CEF,EE,,,,C
@M04916:6:000000000-BRH9H:1:1101:13114:1446 1:N:0:2
TTCCTCAGTTTTACCTACAACACAGAAACAATGATATTACCTACCCCATGGACTGTTGTGAAGATTAAATGAATTAGTACATTTAC
+
8-8AAC--CEEE,CEFCF,@F8B8,,,6C,,C,<E9FF<CCF9FCED,E,,,;E<FF9E,,,,,CE<<,C,,,CF98FAF9EEF<E
";

  import std.stdio;
  import std.algorithm;
  auto tokenizer = input.tokenParser;
  log("First Item");
  auto item1 = tokenizer.nextToken;
  auto id1 = item1.entryid.value(tokenizer.window);
  printFastqEntry(item1.value(tokenizer.window));
  assert( id1 == "M04916:6:000000000-BRH9H:1:1101:14943:1392", "Got: " ~ id1);

  // writefln("Got length: %s, fields: %s" , item1.hfields.length, item1.hfields);
  // writefln("Fields %s", item1.hfields.map!(a => a.value(tokenizer.window)));
  assert(item1.hfields.length == 7);
  assert(item1.hfields[0].value(tokenizer.window) == "M04916");
  assert(item1.hfields[1].value(tokenizer.window) == "6");
  assert(item1.hfields[2].value(tokenizer.window) == "000000000-BRH9H");
  auto extra1 = item1.extras[0].value(tokenizer.window);
  writefln("Extra %s", item1.extras.map!(a => a.value(tokenizer.window)));
  assert( extra1 == "1", "Got: " ~ extra1 ); // "First sub-field of this: 1:N:0:2"
  assert(item1.extras[1].value(tokenizer.window) == "N"); // "Second sub-field of this: 1:N:0:2"
  auto seq = item1.sequence.value(tokenizer.window);
  writefln("Seq1: %s", seq);
  assert(seq[0] == 'T');
  assert(seq[$-1] == 'G');

  auto item2 = tokenizer.nextToken;

  log("Second Item");
  printFastqEntry(item2.value(tokenizer.window));
  auto id2 = item2.entryid.value(tokenizer.window);
  assert( id2 == "M04916:6:000000000-BRH9H:1:1101:13114:1446", "Got: " ~ id2);
  seq = item2.sequence.value(tokenizer.window);
  assert(seq == "TTCCTCAGTTTTACCTACAACACAGAAACAATGATATTACCTACCCCATGGACTGTTGTGAAGATTAAATGAATTAGTACATTTAC");
  auto qual = item2.quality.value(tokenizer.window);
  writeln("Item2 qual: ", qual);
  assert(qual == "8-8AAC--CEEE,CEFCF,@F8B8,,,6C,,C,<E9FF<CCF9FCED,E,,,;E<FF9E,,,,,CE<<,C,,,CF98FAF9EEF<E");
  assert(item2.extras[1].value(tokenizer.window) == "N");
  auto extra2 = item2.extras[0].value(tokenizer.window);
  writefln("Extra2 %s", item2.extras.map!(a => a.value(tokenizer.window)));



  log("Third Item");
  auto item3 = tokenizer.nextToken;
  assert(item3.entryid.length == 0);

  tokenizer.release(item1.endPos);
  item2.release(item1.endPos);
  printFastqEntry(item2);

  auto concrete = item2.value(tokenizer.window);

  assert(concrete.entryid == "M04916:6:000000000-BRH9H:1:1101:13114:1446");
  seq = concrete.sequence;
  assert(seq == "TTCCTCAGTTTTACCTACAACACAGAAACAATGATATTACCTACCCCATGGACTGTTGTGAAGATTAAATGAATTAGTACATTTAC");
}
