import gzip

class Read:
    def __init__(self, name, seq, strand, quality, phred64=False):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.quality = quality
        self.phred64 = phred64


class FastqReader:
    def __init__(self, filename, hasQuality=True, phred64=False):
        self.filename = filename
        self.hasQuality = hasQuality
        self.phred64 = phred64
        self.zipped = False
        self.file = None
        self.buf = bytearray(1024 * 1024)
        self.bufDataLen = 0
        self.bufUsedLen = 0
        self.hasNoLineBreakAtEnd = False
        self.init()

    def init(self):
        if self.filename.endswith(".gz"):
            self.file = gzip.open(self.filename, "rt")
            self.zipped = True
        else:
            self.file = open(self.filename, "rt")
        self.readToBuf()

    def readToBuf(self):
        self.buf = self.file.read(1024 * 1024).encode()
        self.bufDataLen = len(self.buf)
        self.bufUsedLen = 0
        if self.bufDataLen < 1024 * 1024:
            if self.buf[-1] != ord('\n'):
                self.hasNoLineBreakAtEnd = True

    def getLine(self):
        if self.bufUsedLen >= self.bufDataLen and self.eof():
            return None

        start = end = self.bufUsedLen
        while end < self.bufDataLen:
            if self.buf[end] not in b'\r\n':
                end += 1
            else:
                break

        if end < self.bufDataLen or self.bufDataLen < 1024 * 1024:
            line = self.buf[start:end].decode()
            end += 1
            if end < self.bufDataLen - 1 and self.buf[end - 1] == ord('\r') and self.buf[end] == ord('\n'):
                end += 1
            self.bufUsedLen = end
            return line
        else:
            line = self.buf[start:].decode()
            while True:
                self.readToBuf()
                start = 0
                end = 0
                while end < self.bufDataLen:
                    if self.buf[end] not in b'\r\n':
                        end += 1
                    else:
                        break
                if end < self.bufDataLen or self.bufDataLen < 1024 * 1024:
                    line += self.buf[start:end].decode()
                    end += 1
                    if end < self.bufDataLen - 1 and self.buf[end] == ord('\n'):
                        end += 1
                    self.bufUsedLen = end
                    return line
                else:
                    line += self.buf[start:].decode()

    def eof(self):
        return self.bufUsedLen >= self.bufDataLen and self.file.tell() == self.file.seek(0, 2)

    def read(self):
        if self.eof():
            return None

        name = self.getLine()
        while (not name or not name.startswith("@")) and not self.eof():
            name = self.getLine()
        if not name:
            return None

        sequence = self.getLine()
        strand = self.getLine()

        if not self.hasQuality:
            quality = 'K' * len(sequence)
        else:
            quality = self.getLine()
            if len(quality) != len(sequence):
                print(f"ERROR: sequence and quality have different length:\n{name}\n{sequence}\n{strand}\n{quality}")
                return None
        return Read(name, sequence, strand, quality, phred64=self.phred64)

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

    @staticmethod
    def isZipFastq(filename):
        return filename.endswith((".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz"))

    @staticmethod
    def isFastq(filename):
        return filename.endswith((".fastq", ".fq", ".fasta", ".fa"))


class FastqReaderPair:
    def __init__(self, left, right, interleaved=False):
        self.left = left
        self.right = right
        self.interleaved = interleaved

    def read(self):
        left_read = self.left.read()
        if self.interleaved:
            right_read = self.left.read()
        else:
            right_read = self.right.read()

        if left_read is None or right_read is None:
            return None
        else:
            return (left_read, right_read)


def clearLineBreaks(s):
    if s.endswith('\n'):
        s = s[:-1]
    if s.endswith('\r'):
        s = s[:-1]
    return s

# Usage:
# reader = FastqReader("test.fq")
# while True:
#     read = reader.read()
#     if not read:
#         break
#     print(read.name)
#     print(read.seq)
# reader.close()
