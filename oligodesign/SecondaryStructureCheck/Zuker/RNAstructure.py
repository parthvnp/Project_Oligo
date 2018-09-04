"""
RNAstructure.py by Michael Sloma (last update: 3/11/2015)

This is the Python interface to the RNAstructure C++ library
It exports the following classes:
    RNA
    Dynalign_object
    Multilign_object
    Oligowalk_object
    ProbScan

This is mostly a thin wrapper around the C++ classes that has been
automatically generated by SWIG, but it adds some features including
automatic error checking, pretty printing, and some convenience functions.
For object construction,class methods fromString and fromFile have been added
to each class, which are more user-friendly than the C++ constructors.

For further documentation and examples of usage, see the online documentation
at http://rna.urmc.edu
"""

from __future__ import print_function

try:
    from builtins import object
except ImportError:
    from __builtin__ import object

from oligodesign.SecondaryStructureCheck.Zuker import RNAstructure_wrap
from  oligodesign.SecondaryStructureCheck.Zuker import Error_handling as e
import re

#these are named this way to cooperate with python name manglign
filecodes = {"ct":1,"seq":2,"pfs":3,"sav":4}
polymers = {"rna":True,"RNA":True,"Rna":True,"dna":False,"DNA":False,"Dna":False}
invalid_nucs = re.compile(r'[^AUCGTaucgtxX]').search

def sequence_valid(sequence):
    if bool(invalid_nucs(sequence)):
        raise RuntimeError( "unrecognized nucleotides in sequence" )

def polymer_valid(polymer):
    if polymer not in polymers:
        raise RuntimeError( "unknown strand type '%s'" % polymer )

def filetype_valid(filetype):
    if filetype not in filecodes:
        raise RuntimeError( "unknown file type '%s'" % filetype )

def interpret_filetype(filepath):
    """Attempts to automatically interpret the file type of a seq or ct file \
raises an exception if you hand it anything else"""
    buf = open(filepath)
    firstline = buf.readline()
    secondline = buf.readline()
    if re.match(r'\s*[0-9]+\s+\D(\s+[0-9]+){4,4}',secondline):
        return "ct"
    elif re.match(r"^;",firstline.strip()):
        return "seq"
    else:
        raise IOError( "could not automatically determine sequence file type for RNA initialization: \
please provide file type" )

@e.decorate_methods(e.check_for_errors,e.not_excluded)
class RNA(RNAstructure_wrap.RNA):
    """The RNA class provides an entry point for all the single sequence operations of RNAstructure."""
    @e.check_for_errors
    def __init__(self,*args):
        """Initializes an RNA object using a thin wrapper around the C++ constructor (see
        RNAstructure docs for more details.) The class methods RNA.fromFile and
        RNA.fromString provide a better way to initialize RNA objects.
        """
        RNAstructure_wrap.RNA.__init__(self,*args)

    def __iter__(self):
        return self.iterNucs()

    def __len__(self):
        return self.GetSequenceLength()

    def __repr__(self):
        seq = "".join(nuc for nuc in self)
        return "RNA: %s%s" % (seq if len(self)<20 else seq[:15],
                              "" if len(self)<20 else "... %s"%seq[-4:])
    def __str__(self):
        return "".join(nuc for nuc in self)

    @classmethod
    def fromFile(cls,seqfile,filetype=None,backbone="rna"):
        """Create an RNA from an input file.
        Args:
            seqfile(str): A path to the location of the input file
            filetype(str): The format of the input file. Supported options are "seq",
                "ct","sav" (save file from the Fold program), and "pfs" (save file
                from the partition program). If no file type is provided, this function
                will attempt to infer the filetype automatically, and raise an IOError
                if the file cannot be understood.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        if filetype is None:
            filetype = interpret_filetype(seqfile)
        filetype_valid(filetype)
        polymer_valid(backbone)
        return cls(seqfile,filecodes[filetype],polymers[backbone])

    @classmethod
    def fromString(cls,sequence,backbone = "rna"):
        """Create an RNA object from a sequence
        Args:
            sequence(str): The sequence for the RNA. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. If the
                sequence is an RNA, T is treated as U.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        polymer_valid(backbone)
        sequence_valid(sequence)
        return cls(sequence,polymers[backbone])

    def iterNucs(self):
        """iterate over nucleotide identities.
        note that this is has the same behavior as iterating directly
        on the RNA object."""
        return NucIterator(self)

    def iterIndices(self):
        """iterate over 1-indexed nucleotide indices"""
        return range(1,self.GetSequenceLength()+1)

@e.decorate_methods(e.check_for_errors,e.not_excluded)
class HybridRNA(RNAstructure_wrap.HybridRNA):
    @e.check_for_errors
    def __init__(self,*args):
        """Initializes a HybridRNA using a thin wrapper around the C++ constructor (see
        RNAstructure docs for more details.) The class methods HybridRNA.fromFile and
        HybridRNA.fromString provide a better way to initialize HybridRNA objects.
        """
        RNAstructure_wrap.HybridRNA.__init__(self,*args)

    @classmethod
    def fromFile(cls,seqfile1,seqfile2,filetype=None,backbone='rna'):
        """Create a HybridRNA from an input file.
        Args:
            seqfile1(str): A path to the location of the first input file
            seqfile2(str): A path to the location of the second input file
            filetype(str): The format of the input files. Supported options are "seq",
                "ct","sav" (save file from the Fold program), and "pfs" (save file
                from the partition program). If no file type is provided, this function
                will attempt to infer the filetype automatically, and raise an IOError
                if the file cannot be understood.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        if filetype is None:
            filetype1 = interpret_filetype(seqfile1)
            filetype2 = interpret_filetype(seqfile1)
        else:
            filetype_valid(filetype)
            filetype1 = filetype2 = filetype
        polymer_valid(backbone)
        return cls(seqfile1,filecodes[filetype1],seqfile2,
                filecodes[filetype2],polymers[backbone])

    @classmethod
    def fromString(cls,sequence1,sequence2,backbone='rna'):
        """Create a HybridRNA object from sequences
        Args:
            sequence1(str): The sequence of the first RNA. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. If the
                sequence is an RNA, T is treated as U.
            sequence2(str): The sequence of the first RNA. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. If the
                sequence is an RNA, T is treated as U.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        sequence_valid(sequence1)
        sequence_valid(sequence2)
        polymer_valid(backbone)
        return cls(sequence1,sequence2,polymers[backbone])


@e.decorate_methods(e.check_for_errors,e.not_excluded)
class Dynalign_object(RNAstructure_wrap.Dynalign_object):
    """Calculates the optimal structure common to RNAs by simultaneously aligning and folding"""
    @e.check_for_errors
    def __init__(self,*args):
        """Initializes a Dynalign_object using a thin wrapper around the C++ constructor (see
        RNAstructure docs for more details.) The class methods Dynalign_object.fromFile and
        Dynalign_object.fromString provide a better way to initialize RNA objects.
        """
        RNAstructure_wrap.Dynalign_object.__init__(self,*args)

    @classmethod
    def fromFile(cls,seqfile1,seqfile2,filetype=None,backbone='rna'):
        """Create a Dynalign_object from an input file.
        Args:
            seqfile1(str): A path to the location of the first input file
            seqfile2(str): A path to the location of the second input file
            filetype(str): The format of the input files. Supported options are "seq",
                "ct","sav" (save file from the Fold program), and "pfs" (save file
                from the partition program). If no file type is provided, this function
                will attempt to infer the filetype automatically, and raise an IOError
                if the file cannot be understood.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        if filetype is None:
            filetype1 = interpret_filetype(seqfile1)
            filetype2 = interpret_filetype(seqfile1)
        else:
            filetype_valid(filetype)
            filetype1 = filetype2 = filetype
        polymer_valid(backbone)
        return cls(seqfile1,filecodes[filetype1],seqfile2,
                filecodes[filetype2],polymers[backbone])

    @classmethod
    def fromString(cls,sequence1,sequence2,backbone='rna'):
        """Create a Dynalign_object object from sequences
        Args:
            sequence1(str): The sequence of the first RNA. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. If the
                sequence is an RNA, T is treated as U.
            sequence2(str): The sequence of the first RNA. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. If the
                sequence is an RNA, T is treated as U.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        sequence_valid(sequence1)
        sequence_valid(sequence2)
        polymer_valid(backbone)
        return cls(sequence1,sequence2,polymers[backbone])

@e.decorate_methods(e.check_for_errors,e.not_excluded)
class Multilign_object(RNAstructure_wrap.Multilign_object):
    """This class allows calculation of the common structure for two or more sequences
    by repeated pairwise Dynalign calculations."""
    @e.check_for_errors
    def __init__(self,*args):
        """Initializes a Multilign_object using a thin wrapper around the C++ constructor (see
        RNAstructure docs for more details.) The class method Multilign_object.fromFile
        provides a better way to initialize Multilign_objects.
        """
        RNAstructure_wrap.Multilign_object.__init__(self,*args)

    @classmethod
    def fromFile(cls,files,backbone='rna'):
        """Create a Multilign_object from a set of input files
        Args:
            files(str): A list of lists describing the input and output for the calculation.
                Each list should follow the format: [input file, output file, constraints,
                SHAPE]. Each value is a string with a path to the appropriate file. If no
                constraints or SHAPE information are being used, these should be empty
                strings. Input sequences should be in .seq format; output will be in .ct
                format.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        polymer_valid(backbone)
        return cls(files,polymers[backbone],None)

@e.decorate_methods(e.check_for_errors,e.not_excluded)
class Oligowalk_object(RNAstructure_wrap.Oligowalk_object):
    @e.check_for_errors
    def __init__(self,*args):
        """Initializes a Oligowalk_object using a thin wrapper around the C++ constructor (see
        RNAstructure docs for more details.) The class method Oligowalk_object.fromFile
        provides a better way to initialize Oligowalk_objects.
        To create an empty Oligowalk_object for Oligoscreen calculations, call this with no
        arguments
        """
        RNAstructure_wrap.Oligowalk_object.__init__(self,*args)

    @classmethod
    def fromFile(cls,seqfile,filetype=None):
        """Creates an Oligowalk_object from an input file.
        Args:
            seqfile(str): A path to the location of the input file
            filetype(str): The format of the input file. Supported options are "seq",
                "ct","sav" (save file from the Fold program), and "pfs" (save file
                from the partition program). If no file type is provided, this function
                will attempt to infer the filetype automatically, and raise an IOError
                if the file cannot be understood.
        """
        if filetype is None:
            filetype = interpret_filetype(seqfile)
        filetype_valid(filetype)
        return cls(seqfile,filecodes[filetype])

    @classmethod
    def fromString(cls,sequence):
        """Creates an Oligowalk_object from a sequence
        Args:
            sequence(str): The sequence for the Oligowalk_object. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference.
        """
        sequence_valid(sequence)
        return cls(sequence)

@e.decorate_methods(e.check_for_errors,e.not_excluded)
class ProbScan(RNAstructure_wrap.ProbScan):
    """This class is an RNA with additional functionality to calculate loop probabilities
    There are A LOT of methods available because this inherits from RNA. However, the
    following new public methods are specific to the loop probability calculation:
        probability_of_all_hairpins
        probability_of_all_internal_loops
        probability_of_all_helices
        probability_of_hairpin
        probability_of_internal_loop
        probability_of_helix
        probability_of_stack
        probability_of_multibranch_loop
    """

    @e.check_for_errors
    def __init__(self,*args):
        RNAstructure_wrap.ProbScan.__init__(self,*args)

    def __iter__(self):
        return self.iterNucs()

    def __len__(self):
        return self.GetSequenceLength()

    def __repr__(self):
        seq = "".join(nuc for nuc in self)
        return "ProbScan: %s%s" % (seq if len(self)<20 else seq[:15],
                              "" if len(self)<20 else "... %s"%seq[-4:])
    def __str__(self):
        return "".join(nuc for nuc in self)

    @classmethod
    def fromFile(cls,seqfile,filetype=None,backbone="rna"):
        """Create a ProbScan object from an input file.
        Args:
            seqfile(str): A path to the location of the input file
            filetype(str): The format of the input file. Supported options are "seq",
                and "pfs" (save file from the partition program). If no file type is
                provided, this function will attempt to infer the filetype
                automatically, and raise an IOError if the file cannot be understood.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        if filetype is None:
            filetype = interpret_filetype(seqfile)
        filetype_valid(filetype)
        polymer_valid(backbone)
        return cls(seqfile,filetype=='seq',polymers[backbone])

    @classmethod
    def fromString(cls,sequence,backbone = "rna"):
        """Create a ProbScan object from a string containin the sequence
        Args:
            sequence(str): The sequence for the RNA. Input sequence should contain
                A,C,G,T,U,a,c,g,t,u,x,X. Capitalization makes no difference. If the
                sequence is an RNA, T is treated as U.
            backbone(str,optional): The strand type, which determines which set of
                nearest neighbor parameters are used. RNA parameters are used by
                default, unless backbone is 'dna', in which case DNA parameters are used.
        """
        polymer_valid(backbone)
        sequence_valid(sequence)
        return cls(sequence,polymers[backbone])

    def probability_of_all_hairpins(self,threshold=0.01):
        """Finds all possible hairpins with probabilities greater than some threshold
        Args:
            threshold(float,optional): The minimum probability of hairpins to be
                returned. Default is 0.01
        Returns:
            A list of Hairpin objects
        """
        return map(self.Hairpin,RNAstructure_wrap.ProbScan.probability_of_all_hairpins(
            self,3,RNAstructure_wrap.ProbScan.GetSequenceLength(self),threshold))
    def probability_of_all_internal_loops(self,threshold=0.01,mode='both'):
        """Finds all possible internal and/or bulge  loops with probabilities greater than some threshold
        Args:
            threshold(float, optional): The minimum probability of loops to be returned
                Default is 0.01
            mode(str,optional): Describes what type of loops to search for. Allowed values
                are 'internal', 'bulge', and 'both'. Default is 'both'.
        Returns:
            A list of InternalLoop objects
        """
        return map(self.InternalLoop,RNAstructure_wrap.ProbScan.probability_of_all_internal_loops(self,threshold,mode))
    def probability_of_all_helices(self,threshold=0.01,length=1):
        """Finds all possible helices with probabilities greater than some threshold
        Args:
            threshold(float,optional): The minimum probability of helices to be
                returned. Default is 0.01
            length(int,optional): The size of helices to search for, in base pair STACKS
                (a helix with two base pairs has one stack). Default is 1.
        Returns:
            A list of Helix objects
        """
        return map(self.Helix,RNAstructure_wrap.ProbScan.probability_of_all_helices(self,threshold,length))
    def probability_of_multibranch_loop(self,pairs):
        """Calculates the probability of a multibranch loop, provided as a list of pairs
        Args:
            pairs(iterable): A list or tuple of pairs (any iterable with two items)
                that describe the multibranch loop. Example: [(10,40),(12,30),(31,39)]
        Returns:
            The probability of the multibranch loop, as a float
        """
        p = sorted(pairs)
        mbl = RNAstructure_wrap.multibranch_loop(p[0][0],p[0][1])
        for pair in p[1:]:
            i,j = pair
            RNAstructure_wrap.add_branch(mbl,i,j)
        return RNAstructure_wrap.ProbScan.probability_of_multibranch_loop(self,mbl)

    class Hairpin:
        """A representation of a hairpin loop.
        Attributes:
            probability(float): The probability of this loop
            i (int): The position of the 5' nucleotide closing this loop
            j (int): The position of the 3' nucleotide closing this loop
        """
        def __init__(self,hpt):
            self.probability = hpt.probability
            self.i = hpt.i
            self.j = hpt.j
        def __repr__(self):
            return "Hairpin Loop: %0.3f %d %d" % (self.probability,self.i,self.j)
        __string__ = __repr__

    class InternalLoop:
        """A representation of an internal or bulge loop
        Attributes:
            probability(float): The probability of this loop
            i (int): The position of the 5' nucleotide closing this loop on the exterior
            j (int): The position of the 3' nucleotide closing this loop on the exterior
            k (int): The position of the 5' nucleotide closing this loop on the interior
            l (int): The position of the 3' nucleotide closing this loop on the interior
        """
        def __init__(self,hpt):
            self.probability = hpt.probability
            self.i = hpt.i
            self.j = hpt.j
            self.k = hpt.k
            self.l = hpt.l
        def __repr__(self):
            return "Internal Loop: %0.3f %d %d %d %d" % (self.probability,self.i,self.j,self.k,self.l)

    class Helix:
        """A representation of an internal or bulge loop
        Attributes:
            probability(float): The probability of this helix
            i (int): The position of the 5' nucleotide closing this helix on the exterior
            j (int): The position of the 3' nucleotide closing this helix on the exterior
            k (int): The position of the 5' nucleotide closing this helix on the interior
            l (int): The position of the 3' nucleotide closing this helix on the interior
        """
        def __init__(self,hpt):
            self.probability = hpt.probability
            self.i = hpt.i
            self.j = hpt.j
            self.k = hpt.k
            self.l = hpt.l
        def __repr__(self):
            return "Helix: %0.3f %d %d %d %d" % (self.probability,self.i,self.j,self.k,self.l)

    def iterNucs(self):
        """iterate over nucleotide identities.
        note that this is has the same behavior as iterating directly
        on the ProbScan object."""
        return NucIterator(self)

    def iterIndices(self):
        """iterate over 1-indexed nucleotide indices"""
        return range(1,self.GetSequenceLength()+1)

#iterate over nucleotides for RNA and ProbScan
class NucIterator(object):
    def __init__(self,rna):
        self.i = 0
        self.imax = rna.GetSequenceLength()
        self.rna = rna

    def __iter__(self):
        return self

    def __next__(self):
        self.i += 1
        if self.i > self.imax:
            raise StopIteration
        return self.rna.GetNucleotide(self.i)

    next = __next__     # Python2 support
