/**
 * Common utility module that can be used by both FASTA and FASTQ modules
 *
 */

module fastaq.common.utils;

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
