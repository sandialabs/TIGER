
/* 
---------------------------------------------------------------
ARAGORN v1.2.40 Dean Laslett
---------------------------------------------------------------

    ARAGORN (together with ARWEN at last)
    Detects tRNA, mtRNA, and tmRNA genes in nucleotide sequences
    Copyright (C) 2003-2018 Dean Laslett

    A minimum requirement is at least a 32 bit compiler architecture 
    (variable types int and unsigned int are at least 4 bytes long).
    Please report bugs and suggestions for improvements to the authors.

    E-mail: Dean Laslett:  gaiaquark@gmail.com 
            Bjˆrn Canb‰ck: Bjorn.Canback@biol.lu.se

    Version 1.2.40  June 17th, 2019.
    Thanks to Kelly Williams for feedback
    Thanks to Francisco Ossandon for finding many bugs and testing 
    Thanks to Haruo Suzuki for finding bugs
    Thanks to Sascha Steinbiss for fixing bugs


    Please reference the following papers if you use this
    program as part of any published research.

    Laslett, D. and Canback, B. (2004)
    ARAGORN, a program for the detection of transfer RNA and
    transfer-messenger RNA genes in nucleotide sequences.
    Nucleic Acids Research, 32;11-16.

    Laslett, D. and Canback, B. (2008)
    ARWEN: a program to detect tRNA genes in
    metazoan mitochondrial nucleotide sequences.
    Bioinformatics, 24(2); 172-175.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License, (see below).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


               GNU GENERAL PUBLIC LICENSE
               Version 2, June 1991

 Copyright (C) 1989, 1991 Free Software Foundation, Inc.
 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

       Preamble

  The licenses for most software are designed to take away your
freedom to share and change it.  By contrast, the GNU General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  This
General Public License applies to most of the Free Software
Foundation's software and to any other program whose authors commit to
using it.  (Some other Free Software Foundation software is covered by
the GNU Library General Public License instead.)  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
this service if you wish), that you receive source code or can get it
if you want it, that you can change the software or use pieces of it
in new free programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must show them these terms so they know their
rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  Finally, any free program is threatened constantly by software
patents.  We wish to avoid the danger that redistributors of a free
program will individually obtain patent licenses, in effect making the
program proprietary.  To prevent this, we have made it clear that any
patent must be licensed for everyone's free use or not licensed at all.

  The precise terms and conditions for copying, distribution and
modification follow.

   GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License applies to any program or other work which contains
a notice placed by the copyright holder saying it may be distributed
under the terms of this General Public License.  The "Program", below,
refers to any such program or work, and a "work based on the Program"
means either the Program or any derivative work under copyright law:
that is to say, a work containing the Program or a portion of it,
either verbatim or with modifications and/or translated into another
language.  (Hereinafter, translation is included without limitation in
the term "modification".)  Each licensee is addressed as "you".

Activities other than copying, distribution and modification are not
covered by this License; they are outside its scope.  The act of
running the Program is not restricted, and the output from the Program
is covered only if its contents constitute a work based on the
Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does.

  1. You may copy and distribute verbatim copies of the Program's
source code as you receive it, in any medium, provided that you
conspicuously and appropriately publish on each copy an appropriate
copyright notice and disclaimer of warranty; keep intact all the
notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License
along with the Program.

You may charge a fee for the physical act of transferring a copy, and
you may at your option offer warranty protection in exchange for a fee.

  2. You may modify your copy or copies of the Program or any portion
of it, thus forming a work based on the Program, and copy and
distribute such modifications or work under the terms of Section 1
above, provided that you also meet all of these conditions:

    a) You must cause the modified files to carry prominent notices
    stating that you changed the files and the date of any change.

    b) You must cause any work that you distribute or publish, that in
    whole or in part contains or is derived from the Program or any
    part thereof, to be licensed as a whole at no charge to all third
    parties under the terms of this License.

    c) If the modified program normally reads commands interactively
    when run, you must cause it, when started running for such
    interactive use in the most ordinary way, to print or display an
    announcement including an appropriate copyright notice and a
    notice that there is no warranty (or else, saying that you provide
    a warranty) and that users may redistribute the program under
    these conditions, and telling the user how to view a copy of this
    License.  (Exception: if the Program itself is interactive but
    does not normally print such an announcement, your work based on
    the Program is not required to print an announcement.)

These requirements apply to the modified work as a whole.  If
identifiable sections of that work are not derived from the Program,
and can be reasonably considered independent and separate works in
themselves, then this License, and its terms, do not apply to those
sections when you distribute them as separate works.  But when you
distribute the same sections as part of a whole which is a work based
on the Program, the distribution of the whole must be on the terms of
this License, whose permissions for other licensees extend to the
entire whole, and thus to each and every part regardless of who wrote it.

Thus, it is not the intent of this section to claim rights or contest
your rights to work written entirely by you; rather, the intent is to
exercise the right to control the distribution of derivative or
collective works based on the Program.

In addition, mere aggregation of another work not based on the Program
with the Program (or with a work based on the Program) on a volume of
a storage or distribution medium does not bring the other work under
the scope of this License.

  3. You may copy and distribute the Program (or a work based on it,
under Section 2) in object code or executable form under the terms of
Sections 1 and 2 above provided that you also do one of the following:

    a) Accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of Sections
    1 and 2 above on a medium customarily used for software interchange; or,

    b) Accompany it with a written offer, valid for at least three
    years, to give any third party, for a charge no more than your
    cost of physically performing source distribution, a complete
    machine-readable copy of the corresponding source code, to be
    distributed under the terms of Sections 1 and 2 above on a medium
    customarily used for software interchange; or,

    c) Accompany it with the information you received as to the offer
    to distribute corresponding source code.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form with such
    an offer, in accord with Subsection b above.)

The source code for a work means the preferred form of the work for
making modifications to it.  For an executable work, complete source
code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to
control compilation and installation of the executable.  However, as a
special exception, the source code distributed need not include
anything that is normally distributed (in either source or binary
form) with the major components (compiler, kernel, and so on) of the
operating system on which the executable runs, unless that component
itself accompanies the executable.

If distribution of executable or object code is made by offering
access to copy from a designated place, then offering equivalent
access to copy the source code from the same place counts as
distribution of the source code, even though third parties are not
compelled to copy the source along with the object code.

  4. You may not copy, modify, sublicense, or distribute the Program
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense or distribute the Program is
void, and will automatically terminate your rights under this License.
However, parties who have received copies, or rights, from you under
this License will not have their licenses terminated so long as such
parties remain in full compliance.

  5. You are not required to accept this License, since you have not
signed it.  However, nothing else grants you permission to modify or
distribute the Program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and
all its terms and conditions for copying, distributing or modifying
the Program or works based on it.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the
original licensor to copy, distribute or modify the Program subject to
these terms and conditions.  You may not impose any further
restrictions on the recipients' exercise of the rights granted herein.
You are not responsible for enforcing compliance by third parties to
this License.

  7. If, as a consequence of a court judgment or allegation of patent
infringement or for any other reason (not limited to patent issues),
conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot
distribute so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you
may not distribute the Program at all.  For example, if a patent
license would not permit royalty-free redistribution of the Program by
all those who receive copies directly or indirectly through you, then
the only way you could satisfy both it and this License would be to
refrain entirely from distribution of the Program.

If any portion of this section is held invalid or unenforceable under
any particular circumstance, the balance of the section is intended to
apply and the section as a whole is intended to apply in other
circumstances.

It is not the purpose of this section to induce you to infringe any
patents or other property right claims or to contest validity of any
such claims; this section has the sole purpose of protecting the
integrity of the free software distribution system, which is
implemented by public license practices.  Many people have made
generous contributions to the wide range of software distributed
through that system in reliance on consistent application of that
system; it is up to the author/donor to decide if he or she is willing
to distribute software through any other system and a licensee cannot
impose that choice.

This section is intended to make thoroughly clear what is believed to
be a consequence of the rest of this License.

  8. If the distribution and/or use of the Program is restricted in
certain countries either by patents or by copyrighted interfaces, the
original copyright holder who places the Program under this License
may add an explicit geographical distribution limitation excluding
those countries, so that distribution is permitted only in or among
countries not thus excluded.  In such case, this License incorporates
the limitation as if written in the body of this License.

  9. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of this License which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
this License, you may choose any version ever published by the Free Software
Foundation.

  10. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

       NO WARRANTY

  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO
LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

       END OF TERMS AND CONDITIONS


---------------------------------------------------------------
ARAGORN v1.2.40 Dean Laslett
---------------------------------------------------------------
*/


#include <stdio.h>
#include <stdlib.h>

#ifndef SEEK_SET
#define SEEK_SET        0
#define SEEK_CUR        1
#define SEEK_END        2
#endif


#define NOCHAR          '\0'
#define DLIM            '\n'
#define STRLEN          4001
#define STRLENM1        4000
#define SHORTSTRLEN     51
#define SHORTSTRLENM1   50
#define KEYLEN          15
#define NHELPLINE       181  /* was 174 */
#define INACTIVE        2.0e+35
#define IINACTIVE       2000000001L
#define ITHRESHOLD      2000000000L
#define space(c)        (c==' ')||(c=='\t')||(c=='\n')||(c=='\r')
#define sq(pos)         ((pos + d->psmax - 1L) % d->psmax) + 1L
#define itmparam(x,y)   fputc(x,y)



#define FASTA   0
#define GENBANK 1

#define noGENE  -1
#define tRNA    0
#define tmRNA   1
#define srpRNA  2
#define rRNA    3
#define CDS     4
#define NS      6    /* should be one more than number of types of gene */

#define MAXGCMOD   16
#define MAMMAL_MT  2
#define NGENECODE  34    /* previously 26 */
#define METAZOAN_MT      0
#define STANDARD         1
#define VERTEBRATE_MT    2

#define NAMINOACID 27
#define Phe  0
#define Val  1
#define Leu  2
#define Ile  3
#define Cys  4
#define Gly  5
#define Arg  6
#define Ser  7
#define Ala  8
#define Pro  9
#define Thr  10
#define Tyr  11
#define Asp  12
#define His  13
#define Asn  14
#define Met  15
#define Trp  16
#define Glu  17
#define Gln  18
#define Lys  19
#define Stop 20
#define SeC  21
#define Pyl  22

#define INSERT          -2
#define TERM            -1
#define Adenine         0
#define Cytosine        1
#define Guanine         2
#define Thymine         3
#define AMBIG           4
#define NOBASE          5

#define tRNAthresh      132.0
#define mtRNAdtthresh   91.5
#define mtRNAtthresh    83.5
#define mtRNAdthresh    85.0
#define tmRNAthresh     325.0
#define srpRNAthresh    175.0
#define CDSthresh       100.0
#define PSEUDOGENElevel 0.95

#define RIGHT   0
#define UP      1
#define LEFT    2
#define DOWN    3
#define UPRIGHT 4
#define SLANTDR 5
#define SLANTUR 6
#define SLANTUL 7
#define SLANTDL 8
#define SLANT   5

#define MATX 42  
#define MATY 34



#define ASTEM2_EXT       9
#define ASTEM2_EXTD      4                   /* <= ASTEM2_EXT */
#define ASTEM2_EXTE      5                   /* ASTEM2_EXT - ASTEM2_EXTD */
#define MINTSTEM_DIST      (17 + ASTEM2_EXT)
#define MAXTSTEM_DIST      (26 + ASTEM2_EXT)
#define MAXDSTEM_DIST      9
#define MINDSTEM_DIST      8
#define MININTRONLEN    0
#define MAXINTRONLEN    3000
#define MINCTRNALEN     62
#define MAXCTRNALEN     110
#define MINTRNALEN      (MINCTRNALEN + 1)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM2_EXT)
#define MAXETRNALEN     (MAXTRNALEN + MAXINTRONLEN)
#define VARMAX          26 
#define VARMIN          3
#define VARDIFF         23                /* VARMAX - VARMIN */
#define MINTPTSDIST     50
#define MAXTPTSDIST     321
#define TPWINDOW        (MAXTPTSDIST - MINTPTSDIST + 1)
#define MINTPDIST       50
#define MAXTPDIST       250
#define TPDISTWINDOW    (MAXTPDIST - MINTPDIST + 1)
#define MINTAGDIST      12
#define MAXTAGDIST      102
#define TAGWINDOW       MAXTAGDIST - MINTAGDIST
#define MINRNACDIST     (MINTPDIST - 5)
#define MAXRNACDIST     (MAXTPDIST - 5)
#define MAXPPINTRONDIST 250
#define TMPTRAILER      145
#define MINPPASDIST     MINTSTEM_DIST
#define MAXPPASDIST     MAXTSTEM_DIST + MAXPPINTRONDIST
#define MINPPTSTPDIST   MINTSTEM_DIST + MINTPDIST
#define MAXPPTSTPDIST   MAXTSTEM_DIST+ASTEM2_EXT+MAXTPDIST+MAXPPINTRONDIST
#define MAXTMRNALEN     (4 + MAXPPASDIST + MAXTPDIST + MAXTAGDIST + TMPTRAILER)
#define TSWEEP          1000
#define WRAP            2*MAXETRNALEN
#define NPTAG           33
#define MAXAGENELEN     (MAXETRNALEN + MAXTMRNALEN) 

/*
NOTE: If MAXPPINTRONDIST is increased, then validity of MAXTMRNALEN
and MAXETRNALEN must be ensured. WRAP = 2*MAXETRNALEN determines the length
of wseq, which contains the wrap around for circular sequences. This
must remain equal to or more than 2*MAXTMRNALEN and TSWEEP.
*/

#define BASE   0
#define FSTEM  1
#define BSTEM  2

#define NOID   0
#define DLOOP  1
#define DSTEM  2
#define CLOOP  3
#define VAR    4

#define NA          MAXINTRONLEN
#define ND          100
#define NT          200 
#define NH          2000
#define NTH         3000
#define NC          5000 
#define NGFT        5000
#define NTAG        1273
#define NTAGMAX     1300
#define LSEQ        20000
#define ATBOND      2.5
#define mtNA        1500
#define mtND        150 
#define mtNTH       3000 
#define mtNTM       3
#define mtNCDS      200
#define mtNCDSCODON 6000
#define mtGCBOND    0.0
#define mtATBOND    -0.5
#define mtGTBOND    -1.2
#define mtTTBOND    -2.9
#define mtGGBOND    -3.0
#define mtGABOND    -3.0
#define mtNOBOND    -3.0
#define mtBONDSTAB  1.5
#define mtABONDSTAB 2.0
#define mtTSTTSTAB   -2.5
#define mtTERMSTAB  0.01
#define mtSENDSTAB  0.01
#define mtNSTAB     0.1
#define mt3MMSTAB   1.0
#define mtGCPENALTY 0.8
#define mtGCPENALTYD 2.0
#define mt_DRLmaxlength 16
#define mt_TVRLmaxlength 18
#define mtNCLM 3

#define SRRNAMAXLEN 1500
#define SRRNAMINLEN 600
#define LRRNAMINLEN 1200
#define LRRNAMAXLEN 3000

#define srpMAXLEN   650
#define srpUMAXLEN  300
#define srpUMINLEN  100
#define srpDMAXLEN  300
#define srpDMINLEN  100
#define srpNH       200
#define srpNS       500
#define srpMAXHPL   14
#define srpMAXSP    6
#define srpMAXSTEM  6500
#define srpDISPMAX  4*srpMAXLEN
#define srpMAXSPACER 12 
#define srpMAXNISTEMS 10
#define srpNESTMAX  2

#define cdsMAXLEN   3000
#define NCDS        200 
#define NCDSCODON   1000


typedef struct { long start;
                 long stop;
                 int comp;
                 long antistart;
                 long antistop;
                 int genetype;
                 int pseudogene;
                 int permuted;
                 int detected;
                 char species[SHORTSTRLEN]; } annotated_gene;
                 

typedef struct { char filename[80];
                 FILE *f;
                 char seqname[STRLEN];
                 int bugmode;
                 int datatype;
                 double gc;
                 long filepointer;
                 long ps;
                 long psmax;
                 long seqstart;
                 long seqstartoff;
                 long nextseq;
                 long nextseqoff;
                 int ns,nf;
                 long aseqlen; 
                 int nagene[NS];
                 annotated_gene gene[NGFT]; } data_set;


typedef struct { char name[100];
                 int seq[MAXTRNALEN+1];
                 int eseq[MAXETRNALEN+1];
                 int *ps;
                 int nbase;
                 int comp;
                 long start;
                 long stop;
                 int astem1;
                 int astem2;
                 int aatail;
                 int spacer1;
                 int spacer2;
                 int dstem;
                 int dloop;
                 int cstem;
                 int cloop;
                 int intron;
                 int nintron;
                 int anticodon;
                 int var;
                 int varbp;
                 int tstem;
                 int tloop;
                 int genetype;
                 double energy;
                 int asst;
                 int tps;
                 int tpe;
                 int annotation;
                 int annosc;   } gene;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 double energy; } trna_loop;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_loop;

typedef struct { int *pos;
                 int *looppos;
                 int *end;
                 int stem;
                 int loop;
                 int arm;
                 int anticodon;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_cloop;

typedef struct { int *pos;
                 int stem;
                 int loop;
		         int *end;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_tloop;


typedef struct { int *pos;
                 int *end;
                 int stem;
                 int loop;
                 double energy; } trna_dloop;


typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 double energy; } trna_astem;

typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 unsigned int bondtype;
                 double energy; } mt_trna_astem;

typedef struct { int *pos;
                 int comp;
                 int frame;
                 int codon;
                 int win; } mt_cds_codon;

typedef struct { int *pos1;
                 int *pos2;
                 int comp; } mt_cds;

typedef struct { int *pos1;
                 int *pos2;
                 int comp; } mt_rrna;

typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 int loop; } rrna_hairpin;

typedef struct { int *pos1;
                 int *pos2;
                 int stem; } rrna_stem;

typedef struct { int *pos;
                 int comp;
                 int frame;
                 int codon;
                 int win; } cds_codon;

                
typedef struct { char name[50];
                 char tag[50]; } tmrna_tag_entry;


typedef struct { char genetypename[NS][10];
                 FILE *f;
                 int batch;
                 int batchfullspecies;
                 int repeatsn;
                 int trna;
                 int tmrna;
                 int srprna;
                 int cds;
                 int mtrna;
                 int tvloop;
		         int cloop7;
                 int peptide;
                 int geneticcode;
                 int ngcmod;
                 int gcmod[MAXGCMOD];
                 int gcfix;
                 int discrim;
                 int extastem;
                 int tarm;
                 int tagthresh;
                 int tarmlength;
                 int showconfig;
                 int libflag;
                 int verbose;
                 int linear;
                 int both;
                 int reportpseudogenes;
                 int energydisp;
                 int secstructdisp;
                 int seqdisp;
                 int aataildisp;
                 int aataildiv;
                 int sp1max;
		         int sp2min;
		         int sp2max;
		         int mtxdetect;
                 int mtcdsscan;
                 int mtcompov;
		         int matchacceptor;
                 int maxintronlen;
                 int minintronlen;
                 int minintronlenreport;
                 int ioverlay;
		         int ifixedpos;
		         int ireportminintronlen;
                 int tmstrict;
		         int iamismatch;
                 int loffset;
                 int roffset;
                 long start;
                 int comp;
                 int genespace;
                 int srpspace;
                 int ngene[NS];
                 int nps;
                 int annotated;
                 int dispmatch;
                 int updatetmrnatags;
                 int tagend;
                 int trnalenmisthresh;
                 int tmrnalenmisthresh;
                 int nagene[NS];
                 int nafn[NS];
                 int nafp[NS];
                 int natfpd;
                 int natfptv;
                 int lacds;
                 int ldcds;
                 long nabase;
                 double reportpsthresh;
                 double threshlevel;
                 double trnathresh;
                 double ttscanthresh;
                 double ttarmthresh;
                 double tdarmthresh;
                 double tastemthresh;
                 double tascanthresh;
                 double mttthresh;
                 double mtdthresh;
                 double mtdtthresh;
                 double mttarmthresh;
                 double mtdarmthresh;
                 double tmrnathresh;
                 double tmathresh;
                 double tmcthresh;
                 double tmcathresh;
                 double tmrthresh;
                 double srpthresh;
                 double cdsthresh;
                 double eref[NS];
                 int tmrna_struct[200];
               } csw;





/* Basepair matching matrices */

  int lbp[3][6][6] =
   { { { 0,0,1,1,1,0 },
       { 0,0,1,0,1,0 },
       { 1,1,0,1,1,0 },
       { 1,0,1,0,1,0 },
       { 1,1,1,1,1,0 },
       { 0,0,0,0,0,0 } },
     { { 0,0,0,1,1,0 },
       { 0,0,1,0,1,0 },
       { 0,1,0,1,1,0 },
       { 1,0,1,0,1,0 },
       { 1,1,1,1,1,0 },
       { 0,0,0,0,0,0 } },
     { { 0,0,0,1,1,0 },
       { 0,0,1,0,1,0 },
       { 0,1,0,0,1,0 },
       { 1,0,0,0,1,0 },
       { 1,1,1,1,1,0 },
       { 0,0,0,0,0,0 } } };



  int bp[6][6] = { { 0,0,0,1,1,0 },
                   { 0,0,1,0,1,0 },
                   { 0,1,0,1,1,0 },
                   { 1,0,1,0,1,0 },
                   { 1,1,1,1,1,0 },
                   { 0,0,0,0,0,0 } };

  int wbp[6][6] =
   { { 0,0,0,2,2,0 },
     { 0,0,2,0,2,0 },
     { 0,2,0,1,2,0 },
     { 2,0,1,0,2,0 },
     { 2,2,2,2,2,0 },
     { 0,0,0,0,0,0 } };

  int wcbp[6][6] = { { 0,0,0,1,1,0 },
                     { 0,0,1,0,1,0 },
                     { 0,1,0,0,1,0 },
                     { 1,0,0,0,1,0 },
                     { 1,1,1,1,1,0 },
                     { 0,0,0,0,0,0 } };


  int gc[6][6] = { { 0,0,0,0,0,0 },
                   { 0,0,1,0,1,0 },
                   { 0,1,0,0,1,0 },
                   { 0,0,0,0,0,0 },
                   { 0,1,1,0,1,0 },
                   { 0,0,0,0,0,0 } };

  int gt[6][6] = { { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,1,1,0 },
                   { 0,0,1,0,1,0 },
                   { 0,0,1,1,1,0 },
                   { 0,0,0,0,0,0 } };

  int at[6][6] = { { 0,0,0,1,1,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 1,0,0,0,0,0 },
                   { 1,0,0,0,1,0 },
                   { 0,0,0,0,0,0 } };

  int tt[6][6] = { { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,1,1,0 },
                   { 0,0,0,1,1,0 },
                   { 0,0,0,0,0,0 } };

  int stemterm[6][6] = { { 0,0,1,0,1,0 },
                         { 0,0,0,0,0,0 },
                         { 1,0,0,0,1,0 },
                         { 0,0,0,1,1,0 },
                         { 1,0,1,1,1,0 },
                         { 0,0,0,0,0,0 } };

  int aastemterm[6][6] =
   { { 1,0,1,0,1,0 },
     { 0,0,0,0,0,0 },
     { 1,0,0,0,1,0 },
     { 0,0,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int ggstemterm[6][6] =
   { { 0,0,1,0,1,0 },
     { 0,0,0,0,0,0 },
     { 1,0,1,0,1,0 },
     { 0,0,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int assymst[6][6] = { { 0,0,0,0,0,0 },
                        { 0,0,0,0,0,0 },
                        { 1,0,0,0,1,0 },
                        { 0,0,0,1,1,0 },
                        { 1,0,0,1,1,0 },
                        { 0,0,0,0,0,0 } };

  int assymat[6][6] = { { 0,0,0,1,1,0 },
                        { 0,0,0,0,0,0 },
                        { 0,0,0,0,0,0 },
                        { 0,0,0,0,0,0 },
                        { 0,0,0,1,1,0 },
                        { 0,0,0,0,0,0 } };


  int stackbp[6][6] = { { 0,0,0,1,1,0 },
                        { 0,0,1,0,1,0 },
                        { 0,1,0,1,1,0 },
                        { 1,0,1,1,1,0 },
                        { 1,1,1,1,1,0 },
                        { 0,0,0,0,0,0 } };

  int ggstackbp[6][6] =
   { { 0,0,0,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,1,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };


  int ggbp[6][6] =
   { { 0,0,0,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,1,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int gabp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,0,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int assymagbp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,0,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int stembp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int ggstembp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int gastembp[6][6] =
   { { 1,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int vbp[6][6] =
   { { 0,0,1,4,4,0 },
     { 0,0,4,0,4,0 },
     { 1,4,0,2,4,0 },
     { 4,0,2,0,4,0 },
     { 4,4,4,4,4,0 },
     { 0,0,0,0,0,0 } };

  int tandemid[mtNTM][4] =
   { { 3,2,2,3 },
     { 2,3,3,2 },
     { 3,3,3,3 } };

  double tandem_em[mtNTM] = { -0.5,-0.5,2.0 };


  double send_em[6][6] =
   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.5*mtSENDSTAB,0.0,0.5*mtSENDSTAB,0.0 },
     { 0.0,0.5*mtSENDSTAB,0.0,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
     { 0.0,0.5*mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } };


  double ssend_em[6][6] =
   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
     { 0.0,mtSENDSTAB,0.0,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
     { 0.0,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } };



  int neighbour_map[6][6] =
   { { 0,0,1,0,1,0 },
     { 0,0,0,0,0,0 },
     { 1,0,0,0,1,0 },
     { 0,0,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  double neighbour_em[2][6][6] = {
   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } },

   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,mtNSTAB,0.0,mtNSTAB,0.0 },
     { 0.0,mtNSTAB,0.0,0.0,mtNSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,mtNSTAB,mtNSTAB,0.0,mtNSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } } };

  unsigned int btmap[6][6] =
   { { 0x10000,0x10000,0x1000,0x10,0x00000,0x10000 },
     { 0x10000,0x10000,0x1,0x10000,0x00000,0x10000 },
     { 0x1000,0x1,0x10000,0x100,0x00000,0x10000 },
     { 0x10,0x10000,0x100,0x1000,0x00000,0x10000 },
     { 0x00000,0x00000,0x00000,0x00000,0x00000,0x10000 },
     { 0x10000,0x10000,0x10000,0x10000,0x10000,0x10000 } };


  double bem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };



  int mt_discrim[3][64][6] =
   /* metazoan mt */
   {{{ 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 0,0,0,0,0,0 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 0,0,0,0,0,0 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 }},
  /* standard */
    {{ 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 }},
   /* mammal mt */
    {{ 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },
     { 0,0,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
     { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,1,0,0,1,1 },  { 0,0,1,0,1,1 },

     { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,0,1,1 },
     { 0,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,0,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
     { 0,0,0,0,0,0 },  { 1,0,1,1,1,1 },  { 0,0,1,0,1,1 },  { 1,1,1,1,1,1 },

     { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,0,1,1 },
     { 0,1,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,1,0,0,1,1 },  { 0,0,1,0,1,1 },

     { 1,1,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,0,0,0,1,1 },
     { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,0,0,0,0,0 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,1,0,1,1 },  { 1,1,1,1,1,1 },
     { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,1,1,1,1,1 }}};

/* GENETIC CODES (INDEXED BY ANTICODON) */

  char aapolarity[NAMINOACID+1] = "NNNNPNPPNNPNPPPNNPPP***????";
  char aaletter[NAMINOACID+1] = "FVLICGRSAPTYDHNMWEQK***????";
  char aaname[NAMINOACID][20] =
   { "Phe","Val","Leu","Ile","Cys",
     "Gly","Arg","Ser","Ala","Pro",
     "Thr","Tyr","Asp","His","Asn",
     "Met","Trp","Glu","Gln","Lys",
     "Stop",
     "SeC",
     "Pyl",
     "(Arg|Stop|Ser|Gly)",
     "(Ile|Met)",
     "(Stop|Trp)",
     "(Lys|Asn)" };

  char ambig_aaname[4] = "???";

/*
aamap based on NCBI genetic code table (downloaded 26-Apr-2014)
ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
update (downloaded 17-Jun-2019)
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi 
*/

  int aamap[NGENECODE][64] = {
   /* 0. composite metazoan mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,23,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,24,
     25,Gly,Arg,23,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,26 },
   /* 1. standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
  /* 2. vertebrate mt */
  {  Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Stop,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Stop,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 3. yeast mt */
   { Phe,Val,Thr,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Thr,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Thr,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Thr,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 4. mold, protozoan, and coelenterate mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 5. invertebrate mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 6. ciliate */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys },
   /* 7. deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 8. deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 9. echinoderm and flatworm mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Asn },
   /* 10. Euplotid */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Cys,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 11. bacterial and plant chloroplast */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 12. alternate yeast */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Ser,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 13. Ascidian mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Gly,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Gly,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 14. alternate flatworm mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Glu,Gln,Asn },
   /* 15. Blepharisma */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 16. Chlorophycean mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Leu,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 17. deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 18. deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 19. deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 20. deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 21. trematode mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 22. Scenedesmus obliquus mt*/
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Leu,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Stop,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 23. Thraustochytrium mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Ser,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Stop,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 24. Pterobranchia mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Lys,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 25. Gracilibacteria */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Gly,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 26. Pachysolen tannophilus */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Ala,Met,    /* Leu -> Ala */
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 27. Karyorelict */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys,    /* Pyl/Stop -> Gln */
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Arg,    /* SeC/Stop -> Trp */
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys },  /* Stop -> Gln */
   /* 28. Condylostoma */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys,    /* Pyl/Stop -> Gln */
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Arg,    /* SeC/Stop -> Trp */
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys },  /* Stop -> Gln */
   /* 29. Mesodinium */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Tyr,Glu,Gln,Lys,    /* Pyl/Stop -> Tyr */
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Tyr,Glu,Gln,Lys },  /* Stop -> Tyr */
   /* 30. Peritrich */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Glu,Glu,Gln,Lys,    /* Pyl/Stop -> Glu */
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Glu,Glu,Gln,Lys },  /* Stop -> Glu */
   /* 31. Blastocrithidia */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Glu,Glu,Gln,Lys,    /* Pyl/Stop -> Glu */
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Arg,    /* SeC/Stop -> Trp */
     Ser,Ala,Pro,Thr,
     Glu,Glu,Gln,Lys },  /* Stop -> Glu */
   /* 32. vacant -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* 33. Cephalodiscidae Mitochondrial UAA-Tyr */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Lys,    /* Arg -> Lys */
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Ser,    /* SeC/Stop -> Trp, Arg->Ser */
     Ser,Ala,Pro,Thr,
     Tyr,Glu,Gln,Lys }   /* Stop -> Tyr */
 };




/* POINTERS TO DETECTED GENES */

  gene *ts;


/* HELP MENU */

char helpmenu[NHELPLINE][81] =
{
"----------------------------",
"ARAGORN v1.2.40 Dean Laslett",
"----------------------------\n",
"Please reference the following papers if you use this",
"program as part of any published research.\n",
"Laslett, D. and Canback, B. (2004) ARAGORN, a",
"program for the detection of transfer RNA and transfer-messenger",
"RNA genes in nucleotide sequences",
"Nucleic Acids Research, 32;11-16\n",
"Laslett, D. and Canback, B. (2008) ARWEN: a",
"program to detect tRNA genes in metazoan mitochondrial",
"nucleotide sequences",
"Bioinformatics, 24(2); 172-175.\n\n",
"ARAGORN detects tRNA, mtRNA, and tmRNA genes.\n",
"Usage:",
"aragorn -v -e -s -d -c -l -j -a -q -rn -w -ifro<min>,<max> -t -mt -m", 
"        -rp -ps -gc -tv -seq -br -fasta -fo -o <outfile> <filename>\n",
"<filename> is assumed to contain one or more sequences",
"in FASTA or GENBANK format. Results of the search are printed",
"to STDOUT. All switches are optional and case-insensitive.",
"Unless -i is specified, tRNA genes containing introns",
"are not detected.\n",
"    -m            Search for tmRNA genes.",
"    -t            Search for tRNA genes.",
"                  By default, all are detected. If one of",
"                  -m or -t is specified, then the other", 
"                  is not detected unless specified as well.",
"    -mt           Search for Metazoan mitochondrial tRNA genes.",
"                  tRNA genes with introns not detected. -i,-sr switchs",
"                  ignored. Composite Metazoan mitochondrial",
"                  genetic code used.",
"    -mtmam        Search for Mammalian mitochondrial tRNA",
"                  genes. -i switch ignored. -tv switch set.",
"                  Mammalian mitochondrial genetic code used.",
"    -mtx          Same as -mt but low scoring tRNA genes are", 
"                  not reported.",
"    -mtd          Overlapping metazoan mitochondrial tRNA genes", 
"                  on opposite strands are reported.",
"    -gc<num>      Use the GenBank transl_table = <num> genetic code.",
"    -gcstd        Use standard genetic code.",
"    -gcmet        Use composite Metazoan mitochondrial genetic code.",
"    -gcvert       Use Vertebrate mitochondrial genetic code.",
"    -gcinvert     Use Invertebrate mitochondrial genetic code.",
"    -gcyeast      Use Yeast mitochondrial genetic code.",
"    -gcprot       Use Mold/Protozoan/Coelenterate mitochondrial genetic code.",
"    -gcciliate    Use Ciliate genetic code.",
"    -gcflatworm   Use Echinoderm/Flatworm mitochondrial genetic code",
"    -gceuplot     Use Euplotid genetic code.",
"    -gcbact       Use Bacterial/Plant chloroplast genetic code.",
"    -gcaltyeast   Use alternative Yeast genetic code.",
"    -gcascid      Use Ascidian mitochondrial genetic code.",
"    -gcaltflat    Use alternative flatworm mitochondrial genetic code.",
"    -gcblep       Use Blepharisma genetic code.",
"    -gcchloroph   Use Chlorophycean mitochondrial genetic code.",
"    -gctrem       Use Trematode mitochondrial genetic code.",
"    -gcscen       Use Scenedesmus obliquus mitochondrial genetic code.",
"    -gcthraust    Use Thraustochytrium mitochondrial genetic code.",
"    -gcptero      Use Pterobranchia mitochondrial genetic code.",
"    -gcgrac       Use Gracilibacteria genetic code.",
"    -gcpach       Use Pachysolen tannophilus genetic code.",
"    -gckary       Use Karyorelict genetic code.",
"    -gccond       Use Condylostoma genetic code.",
"    -gcmeso       Use Mesodinium genetic code.",
"    -gcperi       Use Peritrich genetic code.",
"    -gcblast      Use Blastocrithidia genetic code.",
"    -gcceph       Use Cephalodiscidae mitochondrial UAA-Tyr genetic code.",
"                  Individual modifications can be appended using",
"    ,BBB=<aa>     B = A,C,G, or T. <aa> is the three letter",
"                  code for an amino-acid. More than one modification",
"                  can be specified. eg -gcvert,aga=Trp,agg=Trp uses",
"                  the Vertebrate mitochondrial code and the codons",
"                  AGA and AGG changed to Tryptophan.",          
"    -c            Assume that each sequence has a circular",
"                  topology. Search wraps around each end.",
"                  Default setting.",
"    -l            Assume that each sequence has a linear",
"                  topology. Search does not wrap.",
"    -d            Double. Search both strands of each",
"                  sequence. Default setting.",
"    -s  or -s+    Single. Do not search the complementary",
"                  (antisense) strand of each sequence.",
"    -sc or -s-    Single complementary. Do not search the sense",
"                  strand of each sequence.",
"    -i            Search for tRNA genes with introns in",
"                  anticodon loop with maximum length 3000",
"                  bases. Minimum intron length is 0 bases.",
"                  Ignored if -m is specified.",
"    -i<max>       Search for tRNA genes with introns in",
"                  anticodon loop with maximum length <max>",
"                  bases. Minimum intron length is 0 bases.",
"                  Ignored if -m is specified.",
"    -i<min>,<max> Search for tRNA genes with introns in",
"                  anticodon loop with maximum length <max>",
"                  bases, and minimum length <min> bases.",
"                  Ignored if -m is specified.",
"    -io           Same as -i, but allow tRNA genes with long",
"                  introns to overlap shorter tRNA genes.",
"    -if           Same as -i, but fix intron between positions",
"                  37 and 38 on C-loop (one base after anticodon).",
"    -ifo          Same as -if and -io combined.",
"    -ir           Same as -i, but report tRNA genes with minimum",
"                  length <min> bases rather than search for", 
"                  tRNA genes with minimum length <min> bases.",
"                  With this switch, <min> acts as an output filter,",
"                  minimum intron length for searching is still 0 bases.",
"    -tv           Do not search for mitochondrial TV replacement",
"                  loop tRNA genes. Only relevant if -mt used.",
"    -c7           Search for tRNA genes with 7 base C-loops only.",
"    -ss           Use the stricter canonical 1-2 bp spacer1 and",
"                  1 bp spacer2. Ignored if -mt set. Default is to",
"                  allow 3 bp spacer1 and 0-2 bp spacer2, which may",
"                  degrade selectivity.",
"    -j            Display 4-base sequence on 3' end of astem",
"                  regardless of predicted amino-acyl acceptor length.",
"    -jr           Allow some divergence of 3' amino-acyl acceptor",
"                  sequence from NCCA.",
"    -jr4          Allow some divergence of 3' amino-acyl acceptor",
"                  sequence from NCCA, and display 4 bases.",
"    -e            Print out score for each reported gene.",
"    -ps           Lower scoring thresholds to 95% of default levels.",
"    -ps<num>      Change scoring thresholds to <num> percent of default levels.",
"    -rp           Flag possible pseudogenes (score < 100 or tRNA anticodon",
"                  loop <> 7 bases long). Note that genes with score < 100",
"                  will not be detected or flagged if scoring thresholds are not",
"                  also changed to below 100% (see -ps switch).",
"    -rp<num>      Flag possible pseudogenes and change score threshold to <num>",
"                  percent of default levels.",
"    -seq          Print out primary sequence.",
"    -br           Show secondary structure of tRNA gene primary sequence,",
"                  or tRNA domain for tmRNA genes, using round brackets.",
"    -svg          Generate SVG image file code for secondary structure.",
"    -fasta        Print out primary sequence in fasta format.",
"    -fo           Print out primary sequence in fasta format only",
"                  (no secondary structure).", 
"    -fon          Same as -fo, with sequence and gene numbering in header.", 
"    -fos          Same as -fo, with no spaces in header.", 
"    -fons         Same as -fo, with sequence and gene numbering, but no spaces.",
"    -v            Verbose. Prints out information during",
"                  search to STDERR.",
"    -a            Print out tRNA domain for tmRNA genes.",
"    -a7           Restrict tRNA astem length to a maximum of 7 bases",
"    -aa           Display message if predicted iso-acceptor species",
"                  does not match species in sequence name (if present).",
"    -amt<num>     Change annotated tRNA length mismatch reporting threshold to",
"                  <num> bases when searching GENBANK files. Default is 10 bases.",
"    -amm<num>     Change annotated tmRNA length mismatch reporting threshold to",
"                  <num> bases when searching GENBANK files. Default is 30 bases.",
"    -q            Dont print configuration line (which switches",
"                  and files were used).",
"    -rn           Repeat sequence name before summary information.",
"    -o <outfile>  Print output to <outfile>. If <outfile>",
"                  already exists, it is overwritten. By default",
"                  all output goes to stdout.",
"    -w            Print out in batch mode.",
"    -wa           Same as -w, but for 6 or 8 base anticodon",
"                  loops, print possible iso-acceptor species",
"                  as ?(<species>|<species>) instead of ???",
"                  For tRNA genes, batch mode output is in the form:\n",
"                  Sequence name",
"                  N genes found",
"                  1 tRNA-<species> [locus 1] <Apos> (nnn)",
"                  i(<intron position>,<intron length>)",
"                            .          ",
"                            .          ",
"                  N tRNA-<species> [Locus N] <Apos> (nnn)",
"                  i(<intron position>,<intron length>)\n",
"                  N is the number of genes found",
"                  <species> is the tRNA iso-acceptor species",
"                  <Apos> is the tRNA anticodon relative position",
"                  (nnn) is the tRNA anticodon base triplet",
"                  i means the tRNA gene has a C-loop intron\n",
"                  For tmRNA genes, output is in the form:\n",
"                  n tmRNA(p) [Locus n] <tag offset>,<tag end offset>",
"                  <tag peptide>\n",
"                  p means the tmRNA gene is permuted",
"    -wunix        Get around problem with some windows gcc compilers",
"                  (found so far in Strawberry Perl and Active Perl)",
"                  when reading Unix files.",
"                  Execution speed may be slower for large files.",
"                  Execution speed will be a lot slower for files",
"                  with many small sequences." 
};



/* tmRNA TAG PEPTIDE DATABASE */
  
tmrna_tag_entry tagdatabase[NTAGMAX] =
   { { "Acaryochloris marina","ANNIVSFARQRTATAVA"},
     { "Accumulibacter phosphatis","ANDERFALAA"},
     { "Acetobacter pasteurianus","ANDNTEVLAVAA"},
     { "Acetobacterium woodii","AKTEKSYGLALAA"},
     { "Acetohalobium arabaticum","ANDNSYALAAA"},
     { "Achromobacter xylosoxidans","ANDERFALAA"},
     { "Acidaminococcus fermentans","ADDSYALAA"},
     { "Acidaminococcus sp. D21","AEDSYALAA"},
     { "Acidimicrobium ferrooxidans","AEPELALAA"},
     { "Acidiphilium cryptum","ANDNFEALAVAA"},
     { "Acidithiobacillus caldus","ANDSNYALAA"},
     { "Acidithiobacillus ferrivorans","ANDSNYALAA"},
     { "Acidithiobacillus ferrooxidans","ANDSNYALAA"},
     { "Acidobacterium capsulatum","ANNNLALAA"},
     { "Acidobacterium Ellin6076","ANTQFAYAA"},
     { "Acidothermus cellulolyticus","ANSSRADFALAA"},
     { "Acidovorax avenae","ANDERFALAA"},
     { "Acidovorax citrulli","ANDERFALAA"},
     { "Acidovorax sp. JS42","ANDERFALAA"},
     { "Acidovorax sp. KKS102","ANDERFALAA"},
     { "Acinetobacter ADP1","ANDETYALAA"},
     { "Acinetobacter baumannii","ANDETYALAA"},
     { "Acinetobacter oleivorans","ANDETYALAA"},
     { "Acinetobacter sp. ADP1","ANDETYALAA"},
     { "Acinetobacter sp. SH024","ANDETYALAA"},
     { "Actinobacillus actinomycetemcomitans","ANDEQYALAA"},
     { "Actinobacillus pleuropneumoniae","ANDEQYALAA"},
     { "Actinobacillus succinogenes","ANDEQYALAA"},
     { "Actinobacillus suis","ANDEQYALAA"},
     { "Actinomyces naeslundii","ADNTRTDFALAA"},
     { "Actinoplanes missouriensis","AKDNSRADFALAA"},
     { "Actinoplanes sp. SE50/110","ANSKFDADQYALAA"},
     { "Actinosynnema mirum","AKSNDQRAFALAA"},
     { "Advenella kashmirensis","ANDESYALAA"},
     { "Aequorivita sublithincola","GENNYALAA"},
     { "Aerococcus urinae","DKNESQSLAFAA"},
     { "Aeromonas hydrophila 1","ANDENYALAA"},
     { "Aeromonas hydrophila 2","ANDENYALAA"},
     { "Aeromonas salmonicida","ANDENYALAA"},
     { "Aeromonas veronii","ANDENYALAA"},
     { "Aggregatibacter actinomycetemcomitans","ANDEQYALAA"},
     { "Aggregatibacter aphrophilus","ANDEQYALAA"},
     { "Agrobacterium fabrum","ANDNNAKEYALAA"},
     { "Agrobacterium radiobacter","ANDNYAEARLAA"},
     { "Agrobacterium sp. H13-3","ANDNNAKEYALAA"},
     { "Agrobacterium tumefaciens 1","ANDNNAKEYALAA"},
     { "Agrobacterium tumefaciens 2","ANDNNAKECALAA"},
     { "Agrobacterium vitis","ANDNNAQGYAVAA"},
     { "Akkermansia muciniphila","AESNDLALAA"},
     { "Alcaligenes faecalis","ANDERFALAA"},
     { "Alcaligenes viscolactis","ANDERFALAA"},
     { "Alcanivorax borkumensis","ANDDSYALAA"},
     { "Alcanivorax dieselolei","ANDDTYALAA"},
     { "Alicycliphilus denitrificans","ANDERFALAA"},
     { "Alicyclobacillus acidocaldarius","GKANRFTTQNKLALAA"},
     { "Aliivibrio salmonicida","ANDENYALAA"},
     { "Alistipes finegoldii","GNNSYALAA"},
     { "Alkalilimnicola ehrlichii","ANDENYALAA"},
     { "Alkaliphilus metalliredigenes","ANDNYSLAAA"},
     { "Alkaliphilus metalliredigens","ANDNYSLAAA"},
     { "Alkaliphilus oremlandii","ANDNYALAA"},
     { "Allochromatium vinosum","ANDDNYALAA"},
     { "alpha proteobacterium","ANESYALAA"},
     { "Alphaproteobacteria SAR-1","ANDELALAA"},
     { "Alteromonas macleodii","ANDETYALAA"},
     { "Alteromonas sp. SN2","ANDENYALAA"},
     { "Aminobacterium colombiense","VNNNNYALAA"},
     { "Ammonifex degensii","ANNERVALAA"},
     { "Amoebophilus asiaticus","GNNQVALAA"},
     { "Amphibacillus xylanus","GKTNNYSLAAA"},
     { "Amycolatopsis mediterranei","ADSSQREFALAA"},
     { "Amycolicicoccus subflavus","ADNAQRSQSDFALAA"},
     { "Anabaena variabilis","ANNIVKFARKDALVAA"},
     { "Anaerobaculum mobile","ANENYALAA"},
     { "Anaerococcus prevotii","ANNNSEANFALAA"},
     { "Anaerolinea thermophila","VRKSGCRSGRSRTERKRAFGP"},
     { "Anaeromyxobacter dehalogenans","ANEPMALAA"},
     { "Anaeromyxobacter sp. Fw109-5","ANEPMALAA"},
     { "Anaeromyxobacter sp. K","ANEPMALAA"},
     { "Anaplasma centrale","ANDDFVAANDNMETAFVAAA"},
     { "Anaplasma marginale","ANDDFVAANDNMETAFVAAA"},
     { "Anaplasma phagocytophilum","ANDDFVAANDNVETAFVAAA"},
     { "Anoxybacillus flavithermus","GKENYALAA"},
     { "Aquifex aeolicus","APEAELALAA"},
     { "Arcanobacterium haemolyticum","ANKQKSDFALAA"},
     { "Arcobacter butzleri","ANNTNYAPAYAKAA"},
     { "Arcobacter nitrofigilis","ANNTNYAPAYAKVA"},
     { "Arcobacter sp. L","ANNTNYAPAYAKAA"},
     { "Aromatoleum aromaticum","ANDERFAVAA"},
     { "Arthrobacter arilaitensis","AESKRTDFALAA"},
     { "Arthrobacter aurescens","AESKRTDFALAA"},
     { "Arthrobacter chlorophenolicus","AESKRTDFALAA"},
     { "Arthrobacter FB24","AKQTRTDFALAA"},
     { "Arthrobacter phenanthrenivorans","AESKRTDFALAA"},
     { "Arthrobacter sp. FB24","AKQTRTDFALAA"},
     { "Arthrobacter sp. Rue61a","AESKRTDFALAA"},
     { "Arthromitus sp. SFB-mouse-Japan","DKNYSLQAA"},
     { "Arthromitus sp. SFB-rat-Yit","DKNYSLQAA"},
     { "Azoarcus BH72","ANDERFALAA"},
     { "Azoarcus EbN1","ANDERFAVAA"},
     { "Azoarcus sp. BH72","ANDERFALAA"},
     { "Azobacteroides pseudotrichonymphae","GENFYALAA"},
     { "Azorhizobium caulinodans","ANDNYAPVAVAA"},
     { "Azospira oryzae","ANDERFAIAA"},
     { "Azospirillum brasilense","ANDNVAPVAVAA"},
     { "Azospirillum lipoferum","ANDNVAQARLAA"},
     { "Azospirillum sp. B510","ANDNVAQARLAA"},
     { "Azotobacter vinelandii","ANDDNYALAA"},
     { "Bacillus amyloliquefaciens","GKTKSFNQNLALAA"},
     { "Bacillus anthracis","GKQNNLSLAA"},
     { "Bacillus atrophaeus","GKTKSFNQNLALAA"},
     { "Bacillus cellulosilyticus","GKQEDNFAFAA"},
     { "Bacillus cereus","GKQNNLSLAA"},
     { "Bacillus clausii","GKENNNFALAA"},
     { "Bacillus coagulans","GKSNTKLALAA"},
     { "Bacillus cytotoxicus","GKQQNNFALAA"},
     { "Bacillus halodurans","GKENNNFALAA"},
     { "Bacillus licheniformis","GKSNQNLALAA"},
     { "Bacillus megaterium","GKSNNNFALAA"},
     { "Bacillus phage","AKLNITNNELQVA"},
     { "Bacillus pumilus","GKTKSFNQNLALAA"},
     { "Bacillus selenitireducens","GKQDNDFALAAA"},
     { "Bacillus stearothermophilus","GKQNYALAA"},
     { "Bacillus subtilis","GKTNSFNQNVALAA"},
     { "Bacillus thuringiensis","GKQNNLSLAA"},
     { "Bacillus weihenstephanensis","GKQNNLSLAA"},
     { "Bacillusphage G","AKLNITNNELQVA"},
     { "Bacteriovorax marinus","AESNFAPAMAA"},
     { "Bacteroides fragilis","GETNYALAA"},
     { "Bacteroides helcogenes","GENNYALAA"},
     { "Bacteroides salanitronis","GNENYALAA"},
     { "Bacteroides thetaiotaomicron","GETNYALAA"},
     { "Bacteroides vulgatus","GNENYALAA"},
     { "Bartonella bacilliformis","ANDNYAEARLAA"},
     { "Bartonella clarridgeiae","ANDNYAEARLIAA"},
     { "Bartonella grahamii","ANDNYAEARLAA"},
     { "Bartonella henselae","ANDNYAEARLAA"},
     { "Bartonella quintana","ANDNYAEARLAA"},
     { "Bartonella tribocorum","ANDNYAEARLAA"},
     { "Baumannia cicadellinicola","ANNSQYESVALAA"},
     { "Bdellovibrio bacteriovorus","GNDYALAA"},
     { "Beijerinckia indica","ANDNYAPVAVAA"},
     { "Belliella baltica","GESNYAMAA"},
     { "Beutenbergia cavernae","ADSKRTDFALAA"},
     { "Bifidobacterium adolescentis","AKSNRTEFALAA"},
     { "Bifidobacterium animalis","AKSNRTEFALAA"},
     { "Bifidobacterium asteroides","AKSNRTEFALAA"},
     { "Bifidobacterium bifidum","AKSNRTEFALAA"},
     { "Bifidobacterium breve","AKSNRTEFALAA"},
     { "Bifidobacterium dentium","AKSNRTEFALAA"},
     { "Bifidobacterium longum","AKSNRTEFALAA"},
     { "Blastococcus saxobsidens","ADSNRADYALAA"},
     { "Blattabacterium sp. (Blaberus giganteus)","GEKEYAFAA"},
     { "Blattabacterium sp. (Blattella germanica) Bge","GEQQYAFAA"},
     { "Blattabacterium sp. (Cryptocercus punctulatus)","GEKQYAFAA"},
     { "Blattabacterium sp. (Mastotermes darwiniensis)","GEKQYAFAA"},
     { "Blattabacterium sp. (Periplaneta americana)","GEKQYAFAA"},
     { "Blochmannia floridanus","AKNKYNEPVALAA"},
     { "Blochmannia pennsylvanicus","ANNTTYRESVALAA"},
     { "Blochmannia vafer","ANYNYNESAALAA"},
     { "Bolidomonas pacifica chloroplast","ANNILAFNRKSLSFA"},
     { "Bordetella avium","ANDERFALAA"},
     { "Bordetella bronchiseptica","ANDERFALAA"},
     { "Bordetella parapertussis","ANDERFALAA"},
     { "Bordetella pertussis","ANDERFALAA"},
     { "Bordetella petrii","ANDERFALAA"},
     { "Borrelia afzelii","AKNNNFTSSNLVMAA"},
     { "Borrelia bissettii","AKNNNFTSSNLVMAA"},
     { "Borrelia burgdorferi","AKNNNFTSSNLVMAA"},
     { "Borrelia crocidurae","AKNNNFTSSDLVMAA"},
     { "Borrelia duttonii","AKNNNFTSSDLVMAA"},
     { "Borrelia garinii","AKNNNFTSSNLVMAA"},
     { "Borrelia hermsii","ARNNNFTSSNLVMAA"},
     { "Borrelia recurrentis","AKNNNFTSSDLVMAA"},
     { "Borrelia turicatae","AKNNNFTSSNLVMAA"},
     { "Brachybacterium faecium","AEPKRTDFALAA"},
     { "Brachyspira hyodysenteriae","ADEYALAA"},
     { "Brachyspira intermedia","ADEYALAA"},
     { "Brachyspira murdochii","ADEYALAA"},
     { "Brachyspira pilosicoli","ADEYALAA"},
     { "Bradyrhizobium japonicum","ANDNFAPVAQAA"},
     { "Bradyrhizobium sp. BTAi1","ANDNFAPVAQAA"},
     { "Bradyrhizobium sp. ORS 278","ANDNFAPVAQAA"},
     { "Bradyrhizobium sp. S23321","ANDNFAPVAQAA"},
     { "Brevibacillus brevis","GNKQLSLAA"},
     { "Brevibacterium linens","AKSNNRTDFALAA"},
     { "Brucella abortus","ANDNNAQGYALAA"},
     { "Brucella canis","ANDNNAQGYALAA"},
     { "Brucella ceti","ANDNNAQGYALAA"},
     { "Brucella melitensis","ANDNNAQGYALAA"},
     { "Brucella ovis","ANDNNAQGYALAA"},
     { "Brucella suis","ANDNNAQGYALAA"},
     { "Buchnera aphidicola 1","ANNKQNYALAA"},
     { "Buchnera aphidicola 2","ANNKQNYALAA"},
     { "Buchnera aphidicola 3","AKQNQYALAA"},
     { "Burkholderia ambifaria","ANDDTFALAA"},
     { "Burkholderia cenocepacia","ANDDTFALAA"},
     { "Burkholderia cepacia","ANDDTFALAA"},
     { "Burkholderia fungorum","ANDDTFALAA"},
     { "Burkholderia gladioli","ANDETFALAA"},
     { "Burkholderia glumae","ANDDTFALAA"},
     { "Burkholderia graminis","ANDDTFALAA"},
     { "Burkholderia mallei","ANDDTFALAA"},
     { "Burkholderia multivorans","ANDDTFALAA"},
     { "Burkholderia phenoliruptrix","ANDDTFALAA"},
     { "Burkholderia phymatum","ANDDTFALAA"},
     { "Burkholderia phytofirmans","ANDETFALAA"},
     { "Burkholderia pseudomallei","ANDDTFALAA"},
     { "Burkholderia rhizoxinica","ANDETYALAA"},
     { "Burkholderia sp. 383","ANDDTFALAA"},
     { "Burkholderia sp. CCGE1001","ANDDTFALAA"},
     { "Burkholderia sp. CCGE1002","ANDDTFALAA"},
     { "Burkholderia sp. YI23","ANDDTFALAA"},
     { "Burkholderia thailandensis","ANDDTFALAA"},
     { "Burkholderia vietnamiensis","ANDDTFALAA"},
     { "Burkholderia xenovorans","ANDDTFALAA"},
     { "Butyrivibrio proteoclasticus","ANDNLALAA"},
     { "Caldicellulosiruptor bescii","ADKAELALAA"},
     { "Caldicellulosiruptor hydrothermalis","ADRTELALAA"},
     { "Caldicellulosiruptor kristjanssonii","ADKAELALAA"},
     { "Caldicellulosiruptor kronotskyensis","ADKAELALAA"},
     { "Caldicellulosiruptor lactoaceticus","ADKAELALAA"},
     { "Caldicellulosiruptor obsidiansis","AEKPQLALAA"},
     { "Caldicellulosiruptor owensensis","AEKPQLALAA"},
     { "Caldicellulosiruptor saccharolyticus","ADKAELALAA"},
     { "Caldilinea aerophila","AKNTGKAFAFGTPATSVALAA"},
     { "Caldisericum exile","ADYSYALAA"},
     { "Calditerrivibrio nitroreducens","ANDEYALAAA"},
     { "Campylobacter coli","ANNVKFAPAYAKAA"},
     { "Campylobacter concisus","ANNVNFAPAYAKAA"},
     { "Campylobacter curvus","ANNVKFAPAYAKAA"},
     { "Campylobacter fetus 2","ANNVKFAPAYAKAA"},
     { "Campylobacter hominis","ANNAKFAPAYAKIA"},
     { "Campylobacter jejuni","ANNVKFAPAYAKAA"},
     { "Campylobacter lari","ANNVKFAPAYAKAA"},
     { "Campylobacter upsaliensis","ANNAKFAPAYAKVA"},
     { "Candidatus atelocyanobacterium thalassa","ANNIVSFKRVAVAA"},
     { "Capnocytophaga canimorsus","GENNYALAA"},
     { "Capnocytophaga ochracea","GENNYALAA"},
     { "Carboxydothermus hydrogenoformans","ANENYALAA"},
     { "Cardinium endosymbiont","VINNSRRCKFVALRKEEEEDDELRMAA"},
     { "Carnobacterium maltaromaticum","AKNNNNSYALAA"},
     { "Carnobacterium sp. 17-4","DKNNNNSYALAA"},
     { "Catenulispora acidiphila","ANKTQLKSQTAYGLAA"},
     { "Catera virion","ATDTDATVTDAEIEAFFAEEAAALV"},
     { "Caulobacter crescentus","ANDNFAEEFAVAA"},
     { "Caulobacter segnis","ANDNFAEEFAVAA"},
     { "Caulobacter sp. K31","ANDNFAEEFAIAA"},
     { "Cellulomonas fimi","ADNKRTDFALAA"},
     { "Cellulomonas flavigena","ADSKRTDFALAA"},
     { "Cellulophaga algicola","GENNYALAA"},
     { "Cellulophaga lytica","GENNYALAA"},
     { "Cellvibrio gilvus","ADSKRTDFALAA"},
     { "Cellvibrio japonicus","ANDDSYALAA"},
     { "Chelativorans sp. BNC1","ANDNYAEARLAA"},
     { "Chitinophaga pinensis","GESNYAMAA"},
     { "Chlamydia muridarum","AEPKAECEIISFADLNDLRVAA"},
     { "Chlamydia psittaci","AEPKAECEIISFSELSEQRLAA"},
     { "Chlamydia trachomatis","AEPKAECEIISFADLEDLRVAA"},
     { "Chlamydophila abortus","AEPKAKCEIISFSELSEQRLAA"},
     { "Chlamydophila caviae","AEPKAECEIISFSDLTEERLAA"},
     { "Chlamydophila felis","AEPKAECEIISFSDLTQERLAA"},
     { "Chlamydophila pecorum","AEPKAECEIISFSDLLVEERVAA"},
     { "Chlamydophila pneumoniae","AEPKAECEIISLFDSVEERLAA"},
     { "Chlamydophila psittaci","AEPKAECEIISFSELSEQRLAA"},
     { "Chloracidobacterium thermophilum","AETQELALAA"},
     { "Chlorobaculum parvum","ADDYSYAMAA"},
     { "Chlorobium chlorochromatii","ADDYSYAMAA"},
     { "Chlorobium limicola","ADDYSYAMAA"},
     { "Chlorobium luteolum","ADDYSYAMAA"},
     { "Chlorobium phaeobacteroides","ADDYSYAMAA"},
     { "Chlorobium phaeovibrioides","ADDYSYAMAA"},
     { "Chlorobium tepidum","ADDYSYAMAA"},
     { "Chloroflexus aggregans","ANNNARVQPRLALAA"},
     { "Chloroflexus aurantiacus","ANTNTRAQARLALAA"},
     { "Chloroherpeton thalassium","ADDYSYAMAA"},
     { "Chromobacterium violaceum","ANDETYALAA"},
     { "Chromohalobacter salexigens","ANDDNYAQGALAA"},
     { "Chroococcidiopsis PCC6712","ANNIVKFERQAVFA"},
     { "Citrobacter koseri","ANDENYALAA"},
     { "Citrobacter rodentium","ANDENYALAA"},
     { "Clavibacter michiganensis","ANNKQSSFVLAA"},
     { "Cloacamonas acidaminovorans","ANNNYALAA"},
     { "Clostridiales genomosp.","ANKNYSYAAA"},
     { "Clostridium acetobutylicum","DNENNLALAA"},
     { "Clostridium acidurici","ANDNYALAA"},
     { "Clostridium beijerinckii","AEDNFALAA"},
     { "Clostridium botulinum","ANDNFALAA"},
     { "Clostridium cellulolyticum","AKNDNFALAAA"},
     { "Clostridium cellulovorans","DENYLLAA"},
     { "Clostridium clariflavum","AENDNYALAAA"},
     { "Clostridium difficile","ADDNFAIAA"},
     { "Clostridium kluyveri","ENDNLALAA"},
     { "Clostridium lentocellum","AEDNLAIAA"},
     { "Clostridium ljungdahlii","ENNNENLALAA"},
     { "Clostridium perfringens","AEDNFALAA"},
     { "Clostridium phytofermentans","ANDNLAYAA"},
     { "Clostridium saccharolyticum","ANNNELALAA"},
     { "Clostridium sp. BNL1100","AKNDNFALAAA"},
     { "Clostridium sp. SY8519","AKEDNFELAMAA"},
     { "Clostridium sticklandii","ANENYALAA"},
     { "Clostridium tetani","ADDNFVLAA"},
     { "Clostridium thermocellum","ANEDNYALAAA"},
     { "Collimonas fungivorans","ANDNSYALAA"},
     { "Colwellia psychrerythraea","ANDDTFALAA"},
     { "Colwellia sp","ANDDTFALAA"},
     { "Comamonas testosteroni","ANDERFALAA"},
     { "Conexibacter woesei","ADSHEYALAA"},
     { "Coprothermobacter proteolyticus","AEPEFALAA"},
     { "Coraliomargarita akajimensis","GEEQFALAA"},
     { "Corallococcus coralloides","ANDNVELALAA"},
     { "Coriobacterium glomerans","GMAQTKIEPTRNPRARRRAQGNRISTGD"},
     { "Corynebacterium aurimucosum","AEKNSQRDYALAA"},
     { "Corynebacterium diphtheriae","AENTQRDYALAA"},
     { "Corynebacterium efficiens","AEKTQRDYALAA"},
     { "Corynebacterium glutamicum","AEKSQRDYALAA"},
     { "Corynebacterium jeikeium","AENTQRDYALAA"},
     { "Corynebacterium kroppenstedtii","AENTQRDYALAA"},
     { "Corynebacterium pseudotuberculosis","AEKTQRDYALAA"},
     { "Corynebacterium resistens","AENTQRDYALAA"},
     { "Corynebacterium ulcerans","AEKTQRDYALAA"},
     { "Corynebacterium urealyticum","AENTQRDYALAA"},
     { "Corynebacterium variabile","AENTQRDYALAA"},
     { "Coxiella burnetii","ANDSNYLQEAYA"},
     { "Croceibacter atlanticus","GENNYALAA"},
     { "Crocosphaera watsonii","ANNIVSFKRVAVAA"},
     { "Cronobacter sakazakii","ANDENYALAA"},
     { "Cronobacter turicensis","ANDENYALAA"},
     { "Cryptobacterium curtum","DNNKSFGRQYALAA"},
     { "Cupriavidus metallidurans","ANDERYALAA"},
     { "Cupriavidus necator","ANDERYALAA"},
     { "Cupriavidus taiwanensis","ANDERYALAA"},
     { "Cyanidioschyzon merolae Chloroplast","ANQILPFSIPVKHLAV"},
     { "Cyanidium caldarium chloroplast","ANNIIEISNIRKPALVV"},
     { "Cyanobium gracile","ANNIVRFSRQAAPVAA"},
     { "Cyanobium sp. PCC 6904","ANNIVRFSRQAAPVAA"},
     { "Cyanobium sp. PCC 7009","ANNIVRFSRQAAPVAA"},
     { "Cyanophora paradoxa chloroplast","ATNIVRFNRKAAFAV"},
     { "Cyanothece sp. ATCC 51142","ANNIVSFKRVAVAA"},
     { "Cyanothece sp. PCC 7424","ANNIVPFARKAAPVAA"},
     { "Cyanothece sp. PCC 7425","ANNIVPFARKAVAVA"},
     { "Cyanothece sp. PCC 7822","ANNIVPFARKSALVAA"},
     { "Cyanothece sp. PCC 8801","ANNIVSFKRVAVAA"},
     { "Cyclobacterium marinum","GESNYAMAA"},
     { "Cycloclasticus sp. P1","ANDDNYAIAA"},
     { "Cytophaga hutchinsonii","GEESYAMAA"},
     { "Dechloromonas agitata","ANDEQFAIAA"},
     { "Dechloromonas aromatica","ANDEQFAIAA"},
     { "Dechlorosoma suillum","ANDERFAIAA"},
     { "Deferribacter desulfuricans","ANDELALAA"},
     { "Dehalococcoides ethenogenes","GERELVLAG"},
     { "Dehalococcoides sp. CBDB1","GERELVLAG"},
     { "Dehalococcoides sp. VS","GERELVLAG"},
     { "Dehalogenimonas lykanthroporepellens","DAKEISAGLERFRRLKLEGREQKAG"},
     { "Deinococcus deserti","GNQNYALAA"},
     { "Deinococcus geothermalis","GNQNYALAA"},
     { "Deinococcus gobiensis","GNQNYALAA"},
     { "Deinococcus maricopensis","GNNNSTTFALAA"},
     { "Deinococcus proteolyticus","GENNYALAA"},
     { "Deinococcus radiodurans","GNQNYALAA"},
     { "Delftia acidovorans","ANDERFALAA"},
     { "Delftia sp. Cs1-4","ANDERFALAA"},
     { "Denitrovibrio acetiphilus","ANNEHTLAAA"},
     { "Desulfarculus baarsii","ADDYNYAVAA"},
     { "Desulfatibacillum alkenivorans","ADDYNYAMAA"},
     { "Desulfitobacterium hafniense","ANDDNYALAA"},
     { "Desulfobacca acetoxidans","ADNYGYALAA"},
     { "Desulfobacterium autotrophicum","ADDYNYAVAA"},
     { "Desulfobacula toluolica","ADDYNYAVAA"},
     { "Desulfobulbus propionicus","ADDYNYALAA"},
     { "Desulfococcus oleovorans","ADDYNYAVAA"},
     { "Desulfohalobium retbaense","ANDYDYALAA"},
     { "Desulfomicrobium baculatum","ANDNYDYAMAA"},
     { "Desulfomonile tiedjei","ANDYEYALAA"},
     { "Desulforudis audaxviator","AKNETYALAA"},
     { "Desulfotalea psychrophila","ADDYNYAVAA"},
     { "Desulfotomaculum acetoxidans","ANNDYALAA"},
     { "Desulfotomaculum carboxydivorans","ANEEYALAA"},
     { "Desulfotomaculum kuznetsovii","ANEEYALAA"},
     { "Desulfotomaculum reducens","ANEEYALAA"},
     { "Desulfotomaculum ruminis","ANEEYALAA"},
     { "Desulfovibrio aespoeensis","ANNDYDYAIAA"},
     { "Desulfovibrio africanus","ANDYNYSLAA"},
     { "Desulfovibrio alaskensis","ANNDYEYAMAA"},
     { "Desulfovibrio desulfuricans","ANNDYDYAYAA"},
     { "Desulfovibrio desulfuricans 2 (G20)","ANNDYEYAMAA"},
     { "Desulfovibrio magneticus","ANDYDYALAA"},
     { "Desulfovibrio salexigens","ANDNYDYAMAA"},
     { "Desulfovibrio vulgaris","ANNYDYALAA"},
     { "Desulfovibrio yellowstonii","ANNELALAA"},
     { "Desulfurispirillum indicum","ANDENVLAAA"},
     { "Desulfurivibrio alkaliphilus","ADDYAYAAAA"},
     { "Desulfurobacterium thermolithotrophum","ANEELALAA"},
     { "Desulfuromonas acetoxidans","ADTDVSYALAA"},
     { "Dichelobacter nodosus","ANDDNYALAA"},
     { "Dickeya dadantii","ANDENFAPAALAA"},
     { "Dickeya zeae","ANDENFAPAALAA"},
     { "Dictyoglomus thermophilum","ANTNLALAA"},
     { "Dictyoglomus turgidum","ANTNLALAA"},
     { "Dinoroseobacter shibae","ANDNRAPVAVAA"},
     { "Dyadobacter fermentans","GESTYAMAA"},
     { "Edwardsiella tarda","ANDENYALAA"},
     { "Eggerthella lenta","GKNNTQSAPALAMAA"},
     { "Eggerthella sp. YY7918","GKNNTQSAPALAMAA"},
     { "Ehrlichia canis","ANDNFVFANDNNSSVAGLVAA"},
     { "Ehrlichia chaffeensis","ANDNFVFANDNNSSANLVAA"},
     { "Ehrlichia ruminantium 1","ANDNFVSANDNNSTANLVAA"},
     { "Ehrlichia ruminantium 2","ANDNFVSANDNNSTANLVAA"},
     { "Elusimicrobium minutum","GNQTELNWATA"},
     { "Emiliania huxleyi chloroplast","ANNILNFNSKLAIA"},
     { "Emticicia oligotrophica","GNTSYAMAA"},
     { "Enterobacter aerogenes","ANDENYALAA"},
     { "Enterobacter cancerogenus","ANDENYALAA"},
     { "Enterobacter cloacae","ANDENYALAA"},
     { "Enterobacter lignolyticus","ANDENYALAA"},
     { "Enterobacter sakazakii","ANDENYALAA"},
     { "Enterobacter sp. 638","ANDENYALAA"},
     { "Enterococcus durans","AKNENNSYALAA"},
     { "Enterococcus faecalis","AKNENNSFALAA"},
     { "Enterococcus faecium","AKNENNSYALAA"},
     { "Enterococcus hirae","AKNENNSYALAA"},
     { "Erwinia amylovora","ANDENFAPAALAA"},
     { "Erwinia billingiae","ANDENYALAA"},
     { "Erwinia carotovora","ANDENYALAA"},
     { "Erwinia chrysanthemi","ANDENFAPAALAA"},
     { "Erwinia pyrifoliae","AKLKYNESVANDGEYELIAAAA"},
     { "Erwinia sp. Ejp617","AKLYNNIPVANDGEFITPALAA"},
     { "Erwinia tasmaniensis","ANDENFAPAALAA"},
     { "Erysipelothrix rhusiopathiae","GNNSLQFAA"},
     { "Erythrobacter litoralis","ANDNEALALAA"},
     { "Escherichia coli","ANDENYALAA"},
     { "Ethanoligenens harbinense","AKDNVIRVNFGRSEEALAA"},
     { "Eubacterium eligens","ANDNLAYAA"},
     { "Eubacterium limosum","AKENRSYGMALAA"},
     { "Eubacterium rectale","AEDNLAYAA"},
     { "Exiguobacterium sibiricum","GKTNTQLAAA"},
     { "Exiguobacterium sp. AT1b","GKTNTQLAAA"},
     { "Ferrimonas balearica","ANDENYALAA"},
     { "Fervidobacterium nodosum","ANEYVPLAA"},
     { "Fervidobacterium pennivorans","ANEYVPLAA"},
     { "Fibrobacter succinogenes","ADENYALAA"},
     { "Filifactor alocis","ANENNLLAA"},
     { "Finegoldia magna","AEDNNFALAA"},
     { "Flavobacteriaceae bacterium","GDQEFALAA"},
     { "Flavobacterium columnare","GENNYALAA"},
     { "Flavobacterium indicum","GENNYALAA"},
     { "Flavobacterium johnsoniae","GENNYALAA"},
     { "Flexibacter litoralis","GESNYAMAA"},
     { "Flexistipes sinusarabici","ANDEFALAAA"},
     { "Fluviicola taffensis","DNTSYALAA"},
     { "Francisella cf.","ANDSNFAAVAKAA"},
     { "Francisella noatunensis","ANDSNFAAVTKAA"},
     { "Francisella novicida","ANDSNFAAVAKAA"},
     { "Francisella philomiragia","ANDSNFAAVAKAA"},
     { "Francisella sp. TX077308","ANDSNFAAVAKAA"},
     { "Francisella tularensis 1","GNKKANRVAANDSNFAAVAKAA"},
     { "Francisella tularensis 2","ANDSNFAAVAKAA"},
     { "Frankia alni","ANKTQPVTPLYALAA"},
     { "Frankia sp. CcI3","ANKTQPTTPTYALAA"},
     { "Frankia sp. EAN1pec","ATKTQPASSTFALAA"},
     { "Frankia sp. EuI1c","ANSEQSATSAYALAA"},
     { "Frankia symbiont","ANKSQSATPRTFALAA"},
     { "Frateuria aurantia","ANDDNYALAA"},
     { "Fremyella diplosiphon","ANNIVKFARKEALVAA"},
     { "Fusobacterium nucleatum 1","GNKDYALAA"},
     { "Fusobacterium nucleatum 2","GNKEYALAA"},
     { "Gallibacterium anatis","ANDENYALAA"},
     { "Gallionella capsiferriformans","ANDENYALAA"},
     { "gamma proteobacterium","ANDESYALAA"},
     { "Gammaproteobacteria SAR-1","ANNYNYSLAA"},
     { "Gardnerella vaginalis","AKSNRTEFALAA"},
     { "Gemmata obscuriglobus","AEPQYSLAA"},
     { "Gemmatimonas aurantiaca","ANNNLALAA"},
     { "Geobacillus kaustophilus","GKQNYALAA"},
     { "Geobacillus sp. WCH70","GKENYALAA"},
     { "Geobacillus sp. Y4.1MC1","GKENYALAA"},
     { "Geobacillus stearothermophilus","GKQNYALAA"},
     { "Geobacillus thermodenitrificans","GKENYALAA"},
     { "Geobacter bemidjiensis","ADNYDYALAA"},
     { "Geobacter daltonii","ADNYDYALAA"},
     { "Geobacter lovleyi","ADNYNTQPVALAA"},
     { "Geobacter metallireducens","ADNYDYAVAA"},
     { "Geobacter sp. M18","ADNYDYALAA"},
     { "Geobacter sp. M21","ADNYDYALAA"},
     { "Geobacter sulfurreducens","ADNYDYAVAA"},
     { "Geobacter uraniireducens","ADNYNYALAA"},
     { "Geodermatophilus obscurus","ADSSQREFALAA"},
     { "Glaciecola nitratireducens","ANDENYALAA"},
     { "Glaciecola sp. 4H-3-7+YE-5","ANDENYALAA"},
     { "Gloeobacter violaceus","ATNNVVPFARARATVAA"},
     { "Gluconacetobacter diazotrophicus","ANDNSEVLAVAA"},
     { "Gluconacetobacter xylinus","ANDNSEVLAVAA"},
     { "Gluconobacter oxydans","ANDNSEVLAVAA"},
     { "Gordonia bronchialis","ADSNQRDYALAA"},
     { "Gordonia polyisoprenivorans","ADKNQRDYALAA"},
     { "Gordonia rubripertincta","ADSNQRDYALAA"},
     { "Gordonia sp. KTR9","ADSNQRDYALAA"},
     { "Gracilaria tenuistipitata chloroplast","AKNNILTLSRRLIYA"},
     { "Gramella forsetii","GENNYALAA"},
     { "Granulibacter bethesdensis","ANDNHEALAVAA"},
     { "Granulicella mallensis","AEPQFALAA"},
     { "Granulicella tundricola","AEPQFALAA"},
     { "Guillardia theta chloroplast","ASNIVSFSSKRLVSFA"},
     { "Haemophilus ducreyi","ANDEQYALAA"},
     { "Haemophilus influenzae","ANDEQYALAA"},
     { "Haemophilus parainfluenzae","ANDEQYALAA"},
     { "Haemophilus parasuis","ANDEQYALAA"},
     { "Haemophilus somnus","ANDEQYALAA"},
     { "Hahella chejuensis","ANDETYALAA"},
     { "Halanaerobium hydrogeniformans","ANDNSYALAAA"},
     { "Halanaerobium praevalens","ANDNNYTLAAA"},
     { "Haliangium ochraceum","ANDNAVALAA"},
     { "Haliscomenobacter hydrossis","GESNYAMAA"},
     { "Halobacillus halophilus","GESNDNLAVAA"},
     { "Halomonas elongata","ANDDNYAQGALAA"},
     { "Halorhodospira halophila","ANDDNYALAA"},
     { "Halothermothrix orenii","ADNNNYALAAA"},
     { "Halothiobacillus neapolitanus","ANDDNYALAA"},
     { "Hamiltonella defensa","AKINKNRPAANGYMPVAALAA"},
     { "Helicobacter acinonychis","VNNTDYAPAYAKVA"},
     { "Helicobacter bizzozeronii","VNNPNYAPNYAKAA"},
     { "Helicobacter cetorum","VNNTNYAPAYAKVA"},
     { "Helicobacter cinaedi","ANNTNYAPVYAKVA"},
     { "Helicobacter felis","VNNPNYAPNYAKAA"},
     { "Helicobacter hepaticus","ANNANYAPAYAKVA"},
     { "Helicobacter mustelae","ANNKNYAPAYAKVA"},
     { "Helicobacter pylori 1","VNNTDYAPAYAKAA"},
     { "Helicobacter pylori 2","VNNTDYAPAYAKAA"},
     { "Helicobacter pylori 3","VNNADYAPAYAKAA"},
     { "Heliobacillus mobilis","AEDNYALAA"},
     { "Heliobacterium modesticaldum","AEENYALAA"},
     { "Herbaspirillum seropedicae","ANDESYALAA"},
     { "Herminiimonas arsenicoxydans","DNSYALAA"},
     { "Herpetosiphon aurantiacus","GKNTFRAPVALAA"},
     { "Hippea maritima","ADTEYALAA"},
     { "Hirschia baltica","ANDNFAEGELLAA"},
     { "Hydrogenophaga palleronii","ANDERFALAA"},
     { "Hyphomicrobium denitrificans","ANDNYAEAALAA"},
     { "Hyphomicrobium sp. MC1","ANDNYAEAALAA"},
     { "Hyphomonas neptunium","ANDNFAEGELLAA"},
     { "Idiomarina loihiensis","ANDDNYALAA"},
     { "Ignavibacterium album","GEYNYALAA"},
     { "Ilyobacter polytropus","ENNNYALAA"},
     { "Intrasporangium calvum","ANSKRTDFALAA"},
     { "Isoptericola variabilis","ADNKRTDFTLAA"},
     { "Jannaschia sp. CCS1","ANDNRAPAMALAA"},
     { "Janthinobacterium sp. Marseille","ANDNSYALAA"},
     { "Jonesia denitrificans","ADTKRTDFALAA"},
     { "Kangiella koreensis","ANEDNYALAA"},
     { "Ketogulonicigenium vulgare","ANNNRAPAMALAA"},
     { "Kineococcus radiotolerans","ADSKRTEFALAA"},
     { "Kitasatospora setae","ANSKRDSQQFALAA"},
     { "Klebsiella oxytoca","ANDENYALAA"},
     { "Klebsiella pneumoniae","ANDENYALAA"},
     { "Kocuria rhizophila","AKSKRTDFALAA"},
     { "Koribacter versatilis","ANTQMAYAA"},
     { "Kosmotoga olearia","ANTEFALAA"},
     { "Kribbella flavida","ADSKRSSFALAA"},
     { "Krokinobacter sp. 4H-3-7-5","GENNYALAA"},
     { "Kyrpidia tusciae","ANKQELALAA"},
     { "Kytococcus sedentarius","ANSKRTDFALAA"},
     { "Lacinutrix sp. 5H-3-7-4","GENNYALAA"},
     { "Lactobacillus acidophilus","ANNKNSYALAA"},
     { "Lactobacillus amylovorus","ANNKNSYALAA"},
     { "Lactobacillus brevis","AKNNNNSYALAA"},
     { "Lactobacillus buchneri","AKNNNNSYALAA"},
     { "Lactobacillus casei","AKNENSYALAA"},
     { "Lactobacillus crispatus","ANNKNSYALAA"},
     { "Lactobacillus delbrueckii 1","AKNENNSYALAA"},
     { "Lactobacillus delbrueckii 2","ANENSYAVAA"},
     { "Lactobacillus fermentum","ANNNSQSYAYAA"},
     { "Lactobacillus gallinarum","ANNKNSYALAA"},
     { "Lactobacillus gasseri","ANNENSYAVAA"},
     { "Lactobacillus helveticus","ANNKNSYALAA"},
     { "Lactobacillus johnsonii","ANNENSYAVAA"},
     { "Lactobacillus kefiranofaciens","ANNKNSYALAA"},
     { "Lactobacillus plantarum","AKNNNNSYALAA"},
     { "Lactobacillus reuteri","ANNNSNSYAYAA"},
     { "Lactobacillus rhamnosus","AKNENSYALAA"},
     { "Lactobacillus ruminis","AKNNNYSYALAA"},
     { "Lactobacillus sakei","ANNNNSYAVAA"},
     { "Lactobacillus salivarius","AKNNNNSYALAA"},
     { "Lactobacillus sanfranciscensis","AKNNNNSYALAA"},
     { "Lactococcus garvieae","AKNNTSYALAA"},
     { "Lactococcus lactis","AKNNTQTYAMAA"},
     { "Lactococcus plantarum","AKNTQTYALAA"},
     { "Lactococcus raffinolactis","AKNTQTYAVAA"},
     { "Laribacter hongkongensis","ANDDTYALAA"},
     { "Lawsonia intracellularis","ANNNYDYALAA"},
     { "Leadbetterella byssophila","GNTSYAMAA"},
     { "Legionella longbeachae","ANDENFAGGEAIAA"},
     { "Legionella pneumophila","ANDENFAGGEAIAA"},
     { "Leifsonia xyli","ANSKSTVSAKADFALAA"},
     { "Leptolyngbya boryana","ANNIVPFARKTAPVAA"},
     { "Leptospira biflexa","ANNEFALAA"},
     { "Leptospira borgpetersenii","ANNELALAA"},
     { "Leptospira interrogans","ANNELALAA"},
     { "Leptospirillum ferriphilum","ANEELALAA"},
     { "Leptospirillum ferrooxidans","ANNEMALAA"},
     { "Leptospirillum groupII","ANEELALAA"},
     { "Leptospirillum groupIII","ANEELALAA"},
     { "Leptospirillum sp. Group II '5-way CG'","ANEELALAA"},
     { "Leptospirillum sp. Group III","ANEELALAA"},
     { "Leptothrix cholodnii","ANDSTYALAA"},
     { "Leptotrichia buccalis","GNDNYALAA"},
     { "Leuconostoc carnosum","AKNENTFAVAA"},
     { "Leuconostoc citreum","AKNENSFAIAA"},
     { "Leuconostoc gasicomitatum","AKNENSFAIAA"},
     { "Leuconostoc gelidum","AKNENSFAIAA"},
     { "Leuconostoc lactis","AKNENSFAIAA"},
     { "Leuconostoc mesenteroides","AKNENSFAIAA"},
     { "Leuconostoc pseudomesenteroides","AKNENSYAIAA"},
     { "Leuconostoc sp. C2","AKNENSFAIAA"},
     { "Liberibacter asiaticus","ANDNSAREVLAA"},
     { "Liberibacter solanacearum","ANDNFAGETRLAA"},
     { "Listeria grayi 1","GKEKQNLAFAA"},
     { "Listeria grayi 2","GKQNNNLAFAA"},
     { "Listeria innocua","GKEKQNLAFAA"},
     { "Listeria ivanovii","GKEKQNLAFAA"},
     { "Listeria monocytogenes","GKEKQNLAFAA"},
     { "Listeria seeligeri","GKEKQNLAFAA"},
     { "Listeria welshimeri","GKEKQNLAFAA"},
     { "Lysinibacillus sphaericus","GKQQNLAFAA"},
     { "Macrococcus caseolyticus","GKTNNFAVAA"},
     { "Magnetococcus marinus","ANDEHYAPAFAAA"},
     { "Magnetococcus sp.","ANDEHYAPAFAAA"},
     { "Magnetospirillum magneticum","ANDNVELAAAA"},
     { "Magnetospirillum magnetotacticum 1","ANDNFAPVAVAA"},
     { "Magnetospirillum magnetotacticum 2","ANDNVELAAAA"},
     { "Mahella australiensis","ADNNAELALAA"},
     { "Mannheimia haemolytica","ANDEQYALAA"},
     { "Mannheimia succiniciproducens","ANDEQYALAA"},
     { "Maribacter sp. HTCC2170","GDNNYALAA"},
     { "Maricaulis maris","ANDNFAEEVALAA"},
     { "Marinithermus hydrothermalis","GNNRYALAA"},
     { "Marinitoga piezophila","AEENYALAA"},
     { "Marinobacter adhaerens","ANDENYALAA"},
     { "Marinobacter aquaeolei","ANDENYALAA"},
     { "Marinobacter hydrocarbonoclasticus","ANDENYALAA"},
     { "Marinobacter sp. BSs20148","ANDENYSLAA"},
     { "Marinomonas mediterranea","ANDENYALAA"},
     { "Marinomonas posidonica","ANDENYALAA"},
     { "Marinomonas sp. MWYL1","ANDENYALAA"},
     { "Marivirga tractuosa","GESNYAMAA"},
     { "Megasphaera elsdenii","AKENNFALAA"},
     { "Meiothermus ruber","GNVRSNSYALAA"},
     { "Meiothermus silvanus","GNTQRSYALAA"},
     { "Melioribacter roseus","GEYNYALAA"},
     { "Melissococcus plutonius","AKKQNYSYAVAA"},
     { "Mesoplasma florum","ANKNEENTNEVPTFMLNAGQANYAFA"},
     { "Mesorhizobium ciceri","ANDNYAEARLAA"},
     { "Mesorhizobium loti","ANDNYAEARLAA"},
     { "Mesorhizobium opportunistum","ANDNYAEARLAA"},
     { "Mesorhizobium sp.","ANDNYAEARLAA"},
     { "Mesostigma viride chloroplast","ANNILPFNRKTAVAV"},
     { "Mesotoga prima","ANNEFALAA"},
     { "Methylacidiphilum infernorum","ANEELALAA"},
     { "Methylibium petroleiphilum","ANDERFALAA"},
     { "Methylobacillus flagellatus","ANDETYALAA"},
     { "Methylobacillus glycogenes","ANDETYALAA"},
     { "Methylobacterium extorquens","ANDNFAPVAVAA"},
     { "Methylobacterium nodulans","ANDNYAPVAVAA"},
     { "Methylobacterium populi","ANDNFAPVAVAA"},
     { "Methylobacterium radiotolerans","ANDNFAPVAVAA"},
     { "Methylobacterium sp. 4-46","ANDNYAPVAVAA"},
     { "Methylocella silvestris","ANDNYAPVAVAA"},
     { "Methylococcus capsulatus","ANDDVYALAA"},
     { "Methylocystis sp. SC2","ANDNYAPVAVAA"},
     { "Methylomicrobium alcaliphilum","ANDENYSMALAA"},
     { "Methylomirabilis oxyfera","ANHELALAA"},
     { "Methylomonas methanica","ANDENYSVALAA"},
     { "Methylophaga sp. JAM1","ANDNNYALAA"},
     { "Methylophaga sp. JAM7","ANDNNYALAA"},
     { "Methylotenera mobilis","ANDETYSLAA"},
     { "Methylotenera versatilis","ANDETYSLAA"},
     { "Methylovorus glucosetrophus","ANDETYALAA"},
     { "Micavibrio aeruginosavorus","ANDNFVVANDNSREAAVAIAA"},
     { "Microbacterium testaceum","ADAKRTDFALAA"},
     { "Microbulbifer degradans","ANDDNYGAQLAA"},
     { "Micrococcus luteus","AESKRTDFALAA"},
     { "Microcystis aeruginosa","ANNIVPFARKAAPVAA"},
     { "Microlunatus phosphovorus","AKSEQRTDFALAA"},
     { "Micromonospora aurantiaca","AKNNRADFALAA"},
     { "Midichloria mitochondrii","ANNKFVPANSDFVPALQAA"},
     { "Mobiluncus curtisii","AERNSTESFALAA"},
     { "Modestobacter marinus","ADSSQRDFALAA"},
     { "Moorella thermoacetica","ADDNLALAA"},
     { "Moranella endobia","ANDSQYESVALAA"},
     { "Moraxella catarrhalis","ANDETYALAA"},
     { "Muricauda ruestringensis","GENNYALAA"},
     { "Mycobacteriophage Bxz1 virion","ATDTDATVTDAEIEAFFAEEAAALV"},
     { "Mycobacterium abscessus","ADSHQRDYALAA"},
     { "Mycobacterium africanum","ADSHQRDYALAA"},
     { "Mycobacterium austroafricanum","ADSNQRDYALAA"},
     { "Mycobacterium avium","ADSHQRDYALAA"},
     { "Mycobacterium bovis","ADSHQRDYALAA"},
     { "Mycobacterium chubuense","ADSNQRDYALAA"},
     { "Mycobacterium gilvum","ADSNQRDYALAA"},
     { "Mycobacterium indicus","ADSHQRDYALAA"},
     { "Mycobacterium intracellulare","ADSHQRDYALAA"},
     { "Mycobacterium leprae","ADSYQRDYALAA"},
     { "Mycobacterium marinum","ADSHQRDYALAA"},
     { "Mycobacterium microti","ADSHQRDYALAA"},
     { "Mycobacterium phage","ATDTDATVTDAEIEAFFAEEAAALV"},
     { "Mycobacterium rhodesiae","ADSNQRDFALAA"},
     { "Mycobacterium smegmatis","ADSNQRDYALAA"},
     { "Mycobacterium sp. MCS","ADTNQRDYALAA"},
     { "Mycobacterium tuberculosis","ADSHQRDYALAA"},
     { "Mycoplasma agalactiae","ANDKKSEEVRVELPAFAIANANANLAFA"},
     { "Mycoplasma arthritidis","GNLETSEDKKLDLQFVMNSQTQQNLLFA"},
     { "Mycoplasma bovis","ANDKKSEEVRLELPAFAIANANANLAFA"},
     { "Mycoplasma capricolum","ANKNEETFEMPAFMMNNASAGANFMFA"},
     { "Mycoplasma conjunctivae","ANKKEDKAVDVNLLASQSFNSNLAFA"},
     { "Mycoplasma crocodyli","GKSKKAENEFSFSNPAFAGNLNLAFA"},
     { "Mycoplasma fermentans","AEDKKAEEVNISSLMIAQKMQSQSNLAFA"},
     { "Mycoplasma gallisepticum","DKTSKELADENFVLNQLASNNYALNF"},
     { "Mycoplasma genitalium 1","DKENNEVLVEPNLIINQQASVNFAFA"},
     { "Mycoplasma genitalium 2","DKENNEVLVDPNLIINQQASVNFAFA"},
     { "Mycoplasma haemofelis","ANKQERESSVVNLLMSQPQDLASLSF"},
     { "Mycoplasma hominis","AEEKQNKQSFVLNQMMSSNPVFAY"},
     { "Mycoplasma hyorhinis","GKENKKEDYSLLMNASTQSNLAFAF"},
     { "Mycoplasma leachii","ANKNEETFEMPAFMMNNASAGANFMFA"},
     { "Mycoplasma mobile","GKEKQLEVSPLLMSSSQSNLVFA"},
     { "Mycoplasma mycoides","ADKNEENFEMPAFMINNASAGANYMFA"},
     { "Mycoplasma penetrans","AKNNKNEAVEVELNDFEINALSQNANLALYA"},
     { "Mycoplasma pneumoniae","DKNNDEVLVDPMLIANQQASINYAFA"},
     { "Mycoplasma pulmonis","GTKKQENDYQDLMISQNLNQNLAFASV"},
     { "Mycoplasma putrefaciens","ANKKTEEFEMPAFMINNASAGANLMFA"},
     { "Mycoplasma synoviae","GNKQSQVEEVTREFSPSLYTFNSNLAYA"},
     { "Myxococcus fulvus","ANDNVELALAA"},
     { "Myxococcus xanthus","ANDNVELALAA"},
     { "Nakamurella multipartita","ADSKRTEFALAA"},
     { "Natranaerobius thermophilus","ADEDYALAAA"},
     { "Nautilia profundicola","AANNTNYSPAVARAAA"},
     { "Neisseria gonorrhoeae","ANDETYALAA"},
     { "Neisseria lactamica","ANDETYALAA"},
     { "Neisseria meningitidis","ANDETYALAA"},
     { "Nephroselmis olivacea chloroplast","TTYHSCLEGHLS"},
     { "Niastella koreensis","GNTQFAMAA"},
     { "Nitratifractor salsuginis","ANNTDYRPAYAHAA"},
     { "Nitratiruptor sp. SB155-2","ANNTDYRPAYAVAA"},
     { "Nitrobacter hamburgensis","ANDNYAPVAQAA"},
     { "Nitrobacter Nb-311A","ANDNYAPVAQAA"},
     { "Nitrobacter winogradskyi","ANDNYAPVAQAA"},
     { "Nitrosococcus halophilus","ANDDNYALAA"},
     { "Nitrosococcus oceani","ANDDNYALAA"},
     { "Nitrosococcus watsonii","ANDDNYALAA"},
     { "Nitrosomonas cryotolerans","ANDENYALAA"},
     { "Nitrosomonas europaea","ANDENYALAA"},
     { "Nitrosomonas eutropha","ANDENYALAA"},
     { "Nitrosomonas sp. AL212","ANDENYALAA"},
     { "Nitrosomonas sp. Is79A3","ANDENYALAA"},
     { "Nitrosospira multiformis","ANDENYALAA"},
     { "Nitrospira defluvii","ANQELALAA"},
     { "Nocardia brasiliensis","ADSNQREYALAA"},
     { "Nocardia cyriacigeorgica","ADSHQREYALAA"},
     { "Nocardia farcinica","ADSHQREYALAA"},
     { "Nocardioides sp. JS614","ANTNRSSFALAA"},
     { "Nocardiopsis alba","ANSKRTEFALAA"},
     { "Nocardiopsis dassonvillei","ANSKRTEFALAA"},
     { "Nostoc azollae","ANNIVKFARREALVAA"},
     { "Nostoc PCC7120","ANNIVKFARKDALVAA"},
     { "Nostoc punctiforme","ANNIVNFARKDALVAA"},
     { "Nostoc sp. PCC 7120","ANNIVKFARKDALVAA"},
     { "Novosphingobium aromaticivorans","ANDNEALALAA"},
     { "Novosphingobium sp. PP1Y","ANDNEALALAA"},
     { "Oceanimonas sp. GK1","ANDENYALAA"},
     { "Oceanithermus profundus","GNDNYALAA"},
     { "Oceanobacillus iheyensis","GKETNQPVLAAA"},
     { "Ochrobactrum anthropi","ANDNKAQGYALAA"},
     { "Odontella sinensis chloroplast","ANNLISSVFKSLSTKQNSLNLSFAV"},
     { "Odoribacter splanchnicus","GENNYALAA"},
     { "Oenococcus oeni","AKNNEPSYALAA"},
     { "Oligotropha carboxidovorans","ANDNYAPVAQAA"},
     { "Olsenella uli","DNDSYQGSYALAA"},
     { "Ornithobacterium rhinotracheale","GNNEYALAA"},
     { "Oscillatoria 6304","ANNIVPFARKAAPVAA"},
     { "Oscillatoria acuminata","ANNIVPFARKAAPVAA"},
     { "Owenweeksia hongkongensis","GENNFALAA"},
     { "Paenibacillus larvae","GKQQNNYALAA"},
     { "Paenibacillus mucilaginosus","GNQKQQLAFAA"},
     { "Paenibacillus polymyxa","GKQQNNYAFAA"},
     { "Paenibacillus sp. JDR-2","GKQQQTYAFAA"},
     { "Paenibacillus sp. Y412MC10","GKQQNNYAFAA"},
     { "Paenibacillus terrae","GKQQNNYAFAA"},
     { "Paludibacter propionicigenes","GENNYALAA"},
     { "Pantoea ananatis","ANDENYALAA"},
     { "Pantoea sp. At-9b","ANDNYYDAPAALAA"},
     { "Pantoea stewartii","ANDENYALAA"},
     { "Pantoea vagans","ANDENYALAA"},
     { "Parabacteroides distasonis","GENNYALAA"},
     { "Parachlamydia acanthamoebae","ADSVSYAAAA"},
     { "Parachlamydia UWE25","ANNSNKIAKVDFQEGTFARAA"},
     { "Paracoccus denitrificans","ANDNRAPVALAA"},
     { "Parvibaculum lavamentivorans","ANDNYAEARLAA"},
     { "Parvularcula bermudensis","ANDNSSEGFALAA"},
     { "Pasteurella multocida","ANDEQYALAA"},
     { "Pavlova lutheri chloroplast","ANNILSFNRVAVA"},
     { "Pectobacterium atrosepticum","ANDENYALAA"},
     { "Pectobacterium carotovora","ANDENYALAA"},
     { "Pectobacterium carotovorum","ANDENYALAA"},
     { "Pectobacterium wasabiae","ANDENYALAA"},
     { "Pediococcus claussenii","AKNNNNSYALAA"},
     { "Pediococcus pentosaceus","AKNNNNSYALAA"},
     { "Pedobacter heparinus","GENNYALAA"},
     { "Pedobacter saltans","ENNYALAA"},
     { "Pelagibacter sp. IMCC9063","ANESYAIAA"},
     { "Pelagibacter ubique","ADESYALAA"},
     { "Pelagibacterium halotolerans","ANDNNKAPVALAA"},
     { "Pelobacter carbinolicus","ADTDVSYALAA"},
     { "Pelobacter propionicus","ADNYNTPVALAA"},
     { "Pelodictyon phaeoclathratiforme","ADDYSYAMAA"},
     { "Pelotomaculum thermopropionicum","AKENYALAA"},
     { "Petrotoga mobilis","GGSSLPKFSWNLA"},
     { "Phaeobacter gallaeciensis","ANDNRAPAMAVAA"},
     { "Photobacterium phosphoreum","ANDENYALAA"},
     { "Photobacterium profundum","ANDENFALAA"},
     { "Photorhabdus asymbiotica","ANDNEYALVA"},
     { "Photorhabdus luminescens","ANDEKYALAA"},
     { "Phycisphaera mikurensis","ANDENTIAGRIGFGNDALRLAA"},
     { "Phytoplasma australiense","GKQTNSASEGDQIYNWVPSQSSQNLQQLAFA"},
     { "Pirellula sp.","AEENFALAA"},
     { "Pirellula staleyi","AESNLALAA"},
     { "Planctomyces brasiliensis","ANKQYAMVA"},
     { "Planctomyces limnophilus","ANTGNYALAA"},
     { "Plectonema boryanum","ANNIVPFARKTAPVAA"},
     { "Polaromonas JS666","ANDERFALAA"},
     { "Polaromonas naphthalenivorans","ANDERFALAA"},
     { "Polaromonas sp. JS666","ANDERFALAA"},
     { "Polymorphum gilvum","ANDNYASDVALAA"},
     { "Polynucleobacter necessarius","ANDERFALAA"},
     { "Porphyra purpurea chloroplast","AENNIIAFSRKLAVA"},
     { "Porphyromonas asaccharolytica","AETRHHPGGRCSEAL"},
     { "Porphyromonas gingivalis","GENNYALAA"},
     { "Prevotella denticola","GENNYALAA"},
     { "Prevotella intermedia","GENNYALAA"},
     { "Prevotella melaninogenica","GENNYALAA"},
     { "Prevotella ruminicola","GNNEYALAA"},
     { "Prochlorococcus marinus 1","ANKIVSFSRQTAPVAA"},
     { "Prochlorococcus marinus 2","ANNIVRFSRQPALVAA"},
     { "Prochlorococcus marinus 3","ANKIVSFSRQTAPVAA"},
     { "Prochlorococcus marinus","ANNIVSFSRQTAPVAA"},
     { "Propionibacterium acidipropionici","ADNKRTDFALAA"},
     { "Propionibacterium acnes 1","AENTRTDFALAA"},
     { "Propionibacterium acnes 2","AENTRTDFALAA"},
     { "Propionibacterium freudenreichii","ADTNRTDFALAA"},
     { "Propionibacterium propionicum","ANNSRTDFALAA"},
     { "Prosthecochloris aestuarii","ADDYSYAMAA"},
     { "Proteobacteria SAR-1, version 1","GENADYALAA"},
     { "Proteobacteria SAR-1, version 2","ANNYNYSLAA"},
     { "Proteobacteria SAR-1, version 3","ADNGYMAAA"},
     { "Proteus mirabilis","ANDNQYKALAA"},
     { "Protochlamydia amoebophila","ANNSNKIAKVDFQEGTFARAA"},
     { "Providencia rettgeri","ANDENYALAA"},
     { "Providencia stuartii","ANDENYALAA"},
     { "Pseudoalteromonas atlantica","ANDENYALAA"},
     { "Pseudoalteromonas haloplanktis","ANDDNYSLAA"},
     { "Pseudoalteromonas sp. SM9913","ANDDNYSLAA"},
     { "Pseudogulbenkiania sp. NH8B","ANDETYALAA"},
     { "Pseudomonas aeruginosa","ANDDNYALAA"},
     { "Pseudomonas brassicacearum","ANDENYGQEFAIAA"},
     { "Pseudomonas chlororaphis","ANDETYGEYALAA"},
     { "Pseudomonas entomophila","ANDENYEGYALAA"},
     { "Pseudomonas fluorescens 1","ANDDQYGAALAA"},
     { "Pseudomonas fluorescens 2","ANDENYGQEFALAA"},
     { "Pseudomonas fluorescens 3 (Pf-5)","ANDETYGDYALAA"},
     { "Pseudomonas fulva","ANDENYEGYALAA"},
     { "Pseudomonas mendocina","ANDDNYALAA"},
     { "Pseudomonas protegens","ANDETYGDYALAA"},
     { "Pseudomonas putida 1","ANDENYGAEYKLAA"},
     { "Pseudomonas stutzeri","ANDDNYEGYALAA"},
     { "Pseudomonas syringae 1","ANDENYGAQLAA"},
     { "Pseudomonas syringae 2","ANDETYGEYALAA"},
     { "Pseudomonas syringae 3","ANDENYGAQLAA"},
     { "Pseudonocardia dioxanivorans","ADKSQRAYALAA"},
     { "Pseudovibrio sp. JE062","ANDNYAMDNAVAA"},
     { "Pseudoxanthomonas spadix","ANDDNYGSDFALAA"},
     { "Pseudoxanthomonas suwonensis","ANDDNYALAA"},
     { "Psychrobacter 2734","ANDENYALAA"},
     { "Psychrobacter arcticus","ANDENYALAA"},
     { "Psychrobacter cryohalolentis","ANDENYALAA"},
     { "Psychrobacter sp. PRwf-1","ANDETYALAA"},
     { "Psychroflexus torquis","GEDNYALAA"},
     { "Psychromonas ingrahamii","ANDSNYSLAA"},
     { "Pusillimonas sp. T7-7","ANDERFALAA"},
     { "Rahnella aquatilis","ANDENYALAA"},
     { "Rahnella sp. Y9602","ANDENYALAA"},
     { "Ralstonia eutropha","ANDERYALAA"},
     { "Ralstonia metallidurans","ANDERYALAA"},
     { "Ralstonia pickettii","ANDERYALAA"},
     { "Ralstonia solanacearum","ANDNRYQLAA"},
     { "Ramlibacter tataouinensis","ANDERFALAA"},
     { "Renibacterium salmoninarum","ANSKRTDFALAA"},
     { "Rhizobium etli","ANDNYAEARLAA"},
     { "Rhizobium leguminosarum","ANDNYAEARLAA"},
     { "Rhodobacter capsulatus","ANDNRAPVALAA"},
     { "Rhodobacter sphaeroides","ANDNRAPVALAA"},
     { "Rhodococcus equi","AESTQREYALAA"},
     { "Rhodococcus erythropolis","ADSNQRDYALAA"},
     { "Rhodococcus jostii","ADSNQRDYALAA"},
     { "Rhodococcus opacus","ADSNQRDYALAA"},
     { "Rhodoferax ferrireducens","ANDERFALAA"},
     { "Rhodomicrobium vannielii","ANDNYAGARPVAIAA"},
     { "Rhodomonas salina","ANNIVPFSRKVALV"},
     { "Rhodopirellula baltica","AEENFALAA"},
     { "Rhodopseudomonas palustris","ANDNYAPVAQAA"},
     { "Rhodopseudomonas palustris 4","ANDNVRMNEVRLAA"},
     { "Rhodospirillum centenum","ANDNTAPALRMAA"},
     { "Rhodospirillum photometricum","ANDNVELAAAA"},
     { "Rhodospirillum rubrum","ANDNVELAAAA"},
     { "Rhodothermus marinus","ANDYSYAMAA"},
     { "Rickettsia africae","ANDNNRSVGHLALAA"},
     { "Rickettsia amblyommii","ANDNNRSVGRLALAA"},
     { "Rickettsia australis","ANDNNRSVDLALAA"},
     { "Rickettsia bellii","ANDNYRSAGTPALAVA"},
     { "Rickettsia conorii","ANDNNRSVGHLALAA"},
     { "Rickettsia heilongjiangensis","ANDNNRSVGRLALAA"},
     { "Rickettsia massiliae","ANDNNRSVGRLALAA"},
     { "Rickettsia montanensis","ANDNNRSVGRLALAA"},
     { "Rickettsia parkeri","ANDNNRSVGHLALAA"},
     { "Rickettsia peacockii","ANDNNRSVGRLALAA"},
     { "Rickettsia philipii","ANDNNRSVGRLALAA"},
     { "Rickettsia prowazekii","ANDNRYVGVPALAAA"},
     { "Rickettsia rhipicephali","ANDNNRSVGRLALAA"},
     { "Rickettsia rickettsii","ANDNNRSVGRLALAA"},
     { "Rickettsia sibirica","ANDNNRSVGHLALAA"},
     { "Rickettsia slovaca","ANDNNRSVGRLALAA"},
     { "Rickettsia typhi","ANDNKRYVGVAALAAA"},
     { "Riemerella anatipestifer","GNEEFALAA"},
     { "Riesia pediculicola","AKTKNYAYAQAA"},
     { "Robiginitalea biformata","GDNNYALAA"},
     { "Roseburia hominis","AEDNLAYAA"},
     { "Roseiflexus castenholzii","ANNNKVVAFKPAMALAA"},
     { "Roseiflexus sp. RS-1","ANTNKVVAFKPAMALAA"},
     { "Roseobacter denitrificans","ANDNRAPVAMAA"},
     { "Roseobacter litoralis","ANDNRAPVAMAA"},
     { "Rothia dentocariosa","AKSKRTDFALAA"},
     { "Rothia mucilaginosa","AESKRTDFALAA"},
     { "Rubrivivax gelatinosus","ANDERFALAA"},
     { "Rubrobacter xylanophilus","ANDREMALAA"},
     { "Ruegeria pomeroyi","ANDNRAPVALAA"},
     { "Ruegeria sp. TM1040","ANDNRAPVALAA"},
     { "Ruminococcus albus","GHGYFAKAS"},
     { "Ruminococcus albus","DNDNFAMAA"},
     { "Runella slithyformis","GEYSYAMAA"},
     { "Ruthia magnifica","ANENNYALAA"},
     { "Saccharomonospora viridis","AKTNSQRDFALAA"},
     { "Saccharophagus degradans","ANDDNYGAQLAA"},
     { "Saccharopolyspora erythraea","ADKSQREFALAA"},
     { "Salinibacter ruber","ADDYSYAMAA"},
     { "Salinispora arenicola","AKQNRADFALAA"},
     { "Salinispora tropica","AKQNRADFALAA"},
     { "Salmonella bongori","ANDENYALAA"},
     { "Salmonella enterica 1","ANDETYALAA"},
     { "Salmonella enterica 2","ANDENYALAA"},
     { "Salmonella enterica 3","ANDETYALAA"},
     { "Salmonella enterica 5","ANDETYALAA"},
     { "Salmonella enterica 6","ANDENYALAA"},
     { "Salmonella paratyphi","ANDENYALAA"},
     { "Salmonella typhimurium","ANDETYALAA"},
     { "Salmonella typhi","ANDETYALAA"},
     { "Sanguibacter keddieii","ADSKRTDFALAA"},
     { "Saprospira grandis","GNTNYALAA"},
     { "Sebaldella termitidis","GNDNYALAA"},
     { "secondary endosymbiont","ANDSQFESKTALAA"},
     { "Segniliparus rotundus","ADTTQRDYALAA"},
     { "Selenomonas ruminantium","DEFDYAYAA"},
     { "Selenomonas sputigena","ANEDYALAA"},
     { "Serratia marcescens","ANDENYALAA"},
     { "Serratia plymuthica","ANDSQFESAALAA"},
     { "Serratia proteamaculans","ANDSQFESAALAA"},
     { "Serratia symbiotica","ANDENYALAA"},
     { "Shewanella amazonensis","ANDDNYALAA"},
     { "Shewanella ANA-3","ANDDNYALAA"},
     { "Shewanella baltica","ANDSNYSLAA"},
     { "Shewanella denitrificans","ANDSNYSLAA"},
     { "Shewanella frigidimarina","ANDSNYSLAA"},
     { "Shewanella halifaxensis","ANDSNYSLAA"},
     { "Shewanella loihica","ANDDNYALAA"},
     { "Shewanella oneidensis","ANDDNYALAA"},
     { "Shewanella pealeana","ANDSNYSLAA"},
     { "Shewanella piezotolerans","ANDDNYSLAA"},
     { "Shewanella putrefaciens","ANDDNYALAA"},
     { "Shewanella PV-4","ANDDNYALAA"},
     { "Shewanella SAR-1","ANDDNYALAA"},
     { "Shewanella SAR-1, version 2","ANNDNYALAA"},
     { "Shewanella SAR-2, version 2","ADYGYMAAA"},
     { "Shewanella sediminis","ANDSNYSLAA"},
     { "Shewanella sp. ANA-3","ANDDNYALAA"},
     { "Shewanella sp. MR-4","ANDDNYALAA"},
     { "Shewanella sp. MR-7","ANDDNYALAA"},
     { "Shewanella sp. W3-18-1","ANDDNYALAA"},
     { "Shewanella violacea","ANDSNYSLAA"},
     { "Shewanella woodyi","ANDDNYALAA"},
     { "Shigella boydii","ANDENYALAA"},
     { "Shigella dysenteriae 1","ANDENYALAA"},
     { "Shigella dysenteriae 2","ANDENYALAA"},
     { "Shigella flexneri","ANDENYALAA"},
     { "Shigella sonnei","ANDENYALAA"},
     { "Shimwellia blattae","ANDENYALAA"},
     { "Sideroxydans lithotrophicus","ANDEKYALAA"},
     { "Silicibacter pomeroyi","ANDNRAPVALAA"},
     { "Silicibacter TM1040","ANDNRAPVALAA"},
     { "Simiduia agarivorans","ANDDNYGAQLAA"},
     { "Simkania negevensis","VDTTEDFYLEAA"},
     { "Sinorhizobium fredii","ANDNYAEARLAA"},
     { "Sinorhizobium medicae","ANDNYAEARLAA"},
     { "Sinorhizobium meliloti","ANDNYAEARLAA"},
     { "Slackia heliotrinireducens","GKSYNTGRMALAA"},
     { "Sodalis glossinidius","ANDSQFESNAALAA"},
     { "Solibacillus silvestris","GKQQNFAFAA"},
     { "Solibacter usitatus","ANTQFAYAA"},
     { "Solitalea canadensis","GENNYALAA"},
     { "Sorangium cellulosum","ANDNAYAVAA"},
     { "Sphaerobacter thermophilus","GNESYALAA"},
     { "Sphaerochaeta coccoides","AKKEDENVSYDAEYAFAA"},
     { "Sphaerochaeta globosa","AKKEDEVSFNAEYAFAA"},
     { "Sphaerochaeta pleomorpha","AKKEDEVSFNAEYALAA"},
     { "Sphingobacterium sp. 21","GENNYALAA"},
     { "Sphingobium chlorophenolicum","ANDNEALALAA"},
     { "Sphingobium japonicum","ANDNEALALAA"},
     { "Sphingobium sp. SYK-6","ANDNEALALAA"},
     { "Sphingomonas elodea","ANDNEALAIAA"},
     { "Sphingomonas wittichii","ANDNEALAIAA"},
     { "Sphingopyxis alaskensis","ANDNEALALAA"},
     { "Spirochaeta africana","AKNEDNVVEVAFGNDDTMLAAA"},
     { "Spirochaeta smaragdinae","ANDADYALAA"},
     { "Spirochaeta thermophila","ANDELALAA"},
     { "Spiroplasma kunkelii","ASKKQKEDKIEMPAFMMNNQLAVSMLAA"},
     { "Spirosoma linguale","GEYNYAMAA"},
     { "Stackebrandtia nassauensis","AKTESRSSFALAA"},
     { "Staphylococcus aureus","GKSNNNFAVAA"},
     { "Staphylococcus carnosus","GKTNNNLAVAA"},
     { "Staphylococcus epidermidis","DKSNNNFAVAA"},
     { "Staphylococcus haemolyticus","DKSNNNFAVAA"},
     { "Staphylococcus lugdunensis","GKSNNNFAVAA"},
     { "Staphylococcus pseudintermedius","GKTNNNFAVAA"},
     { "Staphylococcus saprophyticus","GKENNNFAVAA"},
     { "Staphylococcus xylosus","GKENNNFAVAA"},
     { "Starkeya novella","ANDNYAPVAQAA"},
     { "Stenotrophomonas maltophilia","ANDDNYALAA"},
     { "Stigmatella aurantiaca","DGKDTKANDNVELALAA"},
     { "Streptobacillus moniliformis","GKNNFALAA"},
     { "Streptococcus agalactiae","AKNTNSYALAA"},
     { "Streptococcus bovis","AKNTNSYAVAA"},
     { "Streptococcus constellatus","AKNNNSYALAA"},
     { "Streptococcus criceti","AKNTNSYAVAA"},
     { "Streptococcus dysgalactiae","AKNTNSYALAA"},
     { "Streptococcus equi","AKNNTTYALAA"},
     { "Streptococcus gallolyticus","AKNTNSYAVAA"},
     { "Streptococcus gordonii","AKNNTSYALAA"},
     { "Streptococcus macedonicus","AKNTNSYAVAA"},
     { "Streptococcus mitis","AKNNTSYALAA"},
     { "Streptococcus mutans","AKNTNSYAVAA"},
     { "Streptococcus oralis","AKNNTSYALAA"},
     { "Streptococcus parasanguinis","AKNNNSYALAA"},
     { "Streptococcus parauberis","AKNTNTYALAA"},
     { "Streptococcus pneumoniae","AKNNTSYALAA"},
     { "Streptococcus pseudopneumoniae","AKNNTSYALAA"},
     { "Streptococcus pyogenes","AKNTNSYALAA"},
     { "Streptococcus salivarius","AQLNITAKNTNSYAVAA"},
     { "Streptococcus sanguinis","AKNNNSYALAA"},
     { "Streptococcus sobrinus","AKNTNSYAVAA"},
     { "Streptococcus suis","AKNTNTYALAA"},
     { "Streptococcus thermophilus","AKNTNSYAVAA"},
     { "Streptococcus uberis","AKNTNSYALAA"},
     { "Streptococcus zooepidemicus","AKNNTTYALAA"},
     { "Streptomyces aureofaciens","ANSKRDSQQFALAA"},
     { "Streptomyces avermitilis","ANTKSDSQSFALAA"},
     { "Streptomyces avermitilus","ANTKSDSQSFALAA"},
     { "Streptomyces bingchenggensis","ANTKRDSFALAA"},
     { "Streptomyces cattleya","ANNKRDSFALAA"},
     { "Streptomyces coelicolor","ANTKRDSSQQAFALAA"},
     { "Streptomyces collinus","ANTKRDSSSFALAA"},
     { "Streptomyces flavogriseus","ANSKRDSSAFALAA"},
     { "Streptomyces griseus","ANSKRDSSAFALAA"},
     { "Streptomyces hygroscopicus","ANTKRDSFALAA"},
     { "Streptomyces lividans","ANTKRDSSQQAFALAA"},
     { "Streptomyces scabiei","ANSKSDSPQQQFSLAA"},
     { "Streptomyces sp. SirexAA-E","ANTKRDSSAFALAA"},
     { "Streptomyces thermophilus","AKNTNSYAVAA"},
     { "Streptomyces venezuelae","ANSKSDNSRFALAA"},
     { "Streptomyces violaceusniger","ANTKRDSFALAA"},
     { "Streptosporangium roseum","ANKTHSEVSQGNLALAA"},
     { "Sulcia muelleri","GKKNYALAA"},
     { "Sulfuricurvum kujiense","ANNTNYRPAYAVA"},
     { "Sulfurimonas autotrophica","ANNTNYRPALAVA"},
     { "Sulfurimonas denitrificans","ANNTNYRPAYAVA"},
     { "Sulfurospirillum barnesii","ANNSNYRPAYAVA"},
     { "Sulfurospirillum deleyianum","ANNSNYRPAYALAA"},
     { "Sulfurovum sp. NBC37-1","ANNTDYRPAYAVA"},
     { "Synechococcus elongatus","ANNIVPFARKAAPVAA"},
     { "Synechococcus sp. CC9311","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. CC9605","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. CC9902","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. JA-2-3B'a(2-13)","ANNVVPFARKAAALAA"},
     { "Synechococcus sp. JA-3-3Ab (version 1)","ANNVVPFARKAAALAA"},
     { "Synechococcus sp. JA-3-3Ab (version 2)","ANNVVPFARKAAALAA"},
     { "Synechococcus sp. PCC 6301","ANNIVPFARKAAPVAA"},
     { "Synechococcus sp. PCC 6307","ANNIVRFSRQAAPVAA"},
     { "Synechocystis sp. PCC 6803","ANNIVSFKRVAIAA"},
     { "Synechococcus sp. PCC 6904","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. PCC 7002","ANNIVPFARKAAAVA"},
     { "Synechococcus sp. PCC 7009","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. RCC307","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. WH 7803","ANNIVRFSRQAAPVAA"},
     { "Synechococcus sp. WH 8102","ANNIVRFSRHAAPVAA"},
     { "Syntrophobacter fumaroxidans","ADDYAYAVAA"},
     { "Syntrophomonas wolfei","AEDNFALAA"},
     { "Syntrophothermus lipocalidus","ANNELALAA"},
     { "Syntrophus aciditrophicus","ANDYEYALAA"},
     { "Tannerella forsythensis","GENNYALAA"},
     { "Tannerella forsythia","GENNYALAA"},
     { "Taylorella asinigenitalis","ANDDKFALAA"},
     { "Taylorella equigenitalis","ANDENFALAA"},
     { "Tepidanaerobacter acetatoxydans","ANNDLAYAA"},
     { "Teredinibacter turnerae","ANDDNYGAQLAA"},
     { "Terriglobus roseus","AEPQFALAA"},
     { "Terriglobus saanensis","AEPQFALAA"},
     { "Tetragenococcus halophilus","AKNNNNSYALAA"},
     { "Thalassiosira pseudonana chloroplast","ANNIMPFMFNVVKTNRSLTTLNFAV"},
     { "Thalassiosira weissflogii chloroplast","ANNIIPFIFKAVKTKKEAMALNFAV"},
     { "Thauera sp. MZ1T","ANDERFALAA"},
     { "Thermacetogenium phaeum","ANNEYALAA"},
     { "Thermaerobacter marianensis","ANEELALAA"},
     { "Thermanaerovibrio acidaminovorans","ANDNYALAA"},
     { "Thermincola potens","AEENYALAA"},
     { "Thermoanaerobacter italicus","ADRELAYAA"},
     { "Thermoanaerobacter mathranii","ADRELAYAA"},
     { "Thermoanaerobacter pseudethanolicus","ADRELAYAA"},
     { "Thermoanaerobacter sp. X514","ADRELAYAA"},
     { "Thermoanaerobacter tengcongensis","ADRELAYAA"},
     { "Thermoanaerobacter wiegelii","ADRELAYAA"},
     { "Thermoanaerobacterium saccharolyticum","ANDNLAYAA"},
     { "Thermoanaerobacterium thermosaccharolyticum","ANNDNLAYAA"},
     { "Thermoanaerobacterium xylanolyticum","ANDNLAYAA"},
     { "Thermobaculum terrenum","ANTEYALAA"},
     { "Thermobifida fusca","ANSKRTEFALAA"},
     { "Thermobispora bispora","ANKKHAEVSQASLALAA"},
     { "Thermodesulfatator indicus","ADEYNYAMAA"},
     { "Thermodesulfobacterium commune","ANEYAYALAA"},
     { "Thermodesulfobacterium geofontis","ADEYSYALAA"},
     { "Thermodesulfobium narugense","ANNNSLALAA"},
     { "Thermodesulfovibrio yellowstonii","ANNELALAA"},
     { "Thermomicrobium roseum","GERELALAA"},
     { "Thermomonospora curvata","ANKKQSEFALAA"},
     { "Thermosediminibacter oceani","ANEELALAA"},
     { "Thermosipho africanus","ANEELALAA"},
     { "Thermosipho melanesiensis","ANEEIALAA"},
     { "Thermosynechococcus elongatus","ANNIVPFARKAAAVA"},
     { "Thermotoga lettingae","ANNELALAA"},
     { "Thermotoga maritima ","ANEPVAVAA"},
     { "Thermotoga neapolitana","ANEPVAVAA"},
     { "Thermotoga petrophila","ANEPVAVAA"},
     { "Thermotoga sp. RQ2","ANEPVAVAA"},
     { "Thermotoga thermarum","ANEELALAA"},
     { "Thermovibrio ammonificans","ADETLALAA"},
     { "Thermovirga lienii","ANENYALAA"},
     { "Thermus oshimai","ANKPAYALAA"},
     { "Thermus scotoductus","ANKPAYALAA"},
     { "Thermus sp. CCB_US3_UF1","ANKPAYALAA"},
     { "Thermus thermophilus","ANTNYALAA"},
     { "Thioalkalimicrobium cyclicum","ANDDNYALAA"},
     { "Thioalkalivibrio sp. K90mix","ANDDNYALAA"},
     { "Thiobacillus denitrificans","AKSKAARRNPACSAGVMELKA"},
     { "Thiocystis violascens","ANDDNYALAA"},
     { "Thiomicrospira crunogena","ANDDNYALAA"},
     { "Thiomonas intermedia","ANDSSYALAA"},
     { "Thiomonas sp. 3As","ANDSSYALAA"},
     { "Tistrella mobilis","ANDNRVALAA"},
     { "Tolumonas auensis","ANDETYALAA"},
     { "Tremblaya princeps 1 (Dysmicoccus)","APSNRFTIVANDCIDALVRRAVV"},
     { "Treponema azotonutricium","ADNDNYNYALAA"},
     { "Treponema brennaborense","AEDNRQFALAA"},
     { "Treponema caldaria","ADNDSYALAA"},
     { "Treponema denticola","AENNDSFDYALAA"},
     { "Treponema pallidum","ANSDSFDYALAA"},
     { "Treponema primitia","ANNDSYAFAA"},
     { "Treponema succinifaciens","AKRREDEQSENEQFALAA"},
     { "Trichodesmium erythraeum","ANNIVPFARKQVAALA"},
     { "Tropheryma whipplei","ANLKRTDLSLAA"},
     { "Truepera radiovictrix","GNSNSYALAA"},
     { "Tsukamurella paurometabola","ADSNQRDFALAA"},
     { "Turneriella parva","AENETYALAA"},
     { "uncultured bacterium","ANDNFAPVAVAA"},
     { "Uncultured ciona","ANDEFFDARLRA"},
     { "Uncultured FS1","ANDETYALAA"},
     { "Uncultured FS2","ANDENYALAA"},
     { "Uncultured LEM1","ANDETYALAA"},
     { "Uncultured LEM2","ANDETHALAA"},
     { "Uncultured marineEBAC20E09","ANNDNYALAA"},
     { "Uncultured phakopsora","ANDNSYALAA"},
     { "Uncultured QL1","ANVENYALAA"},
     { "Uncultured RCA1","ANDENYALAA"},
     { "Uncultured RCA2","SNDENYALAA"},
     { "Uncultured RCA4","ANDETYALAA"},
     { "Uncultured remanei","ANDESYALAA"},
     { "Uncultured stronglyoides1","ANDERFALAA"},
     { "Uncultured U01a","ANDSNYALAA"},
     { "Uncultured U02","ANDEQFALAA"},
     { "Uncultured U04","ANDETYALAA"},
     { "Uncultured VLS13","ANDENYALAA"},
     { "Uncultured VLS1","ANDENYALAA"},
     { "Uncultured VLS5","ANDETYALAA"},
     { "Uncultured VLS6","ANDENYALAA"},
     { "Uncultured VLS7","ANDENYALAA"},
     { "Uncultured VLS9","ANDENYALAA"},
     { "Uncultured VLW1","ANDENYALAA"},
     { "Uncultured VLW2","ANDENYALAA"},
     { "Uncultured VLW3","ANDENYALAA"},
     { "Uncultured VLW5","ANDENYALAA"},
     { "Uncultured WW10","ANDENYALAV"},
     { "Uncultured WW11","ANDDNYALAA"},
     { "Uncultured WW1","ANDENYALAA"},
     { "Uncultured WW2","ANDENYALAA"},
     { "Uncultured WW4","ANDGNYALAA"},
     { "Uncultured WW5","ANDENYALAA"},
     { "Uncultured WW7","ANDENCALAA"},
     { "Uncultured WW8","ANDENYALAA"},
     { "Uncultured WW9","ANDENYALAA"},
     { "Ureaplasma parvum","AENKKSSEVELNPAFMASATNANYAFAY"},
     { "Ureaplasma urealyticum","AENKKSSEVELNPAFMASATNANYAFAY"},
     { "Variovorax paradoxus","ANDERFALAA"},
     { "Veillonella parvula","AEENFALAA"},
     { "Verminephrobacter eiseniae","ANDERFALAA"},
     { "Verrucomicrobium spinosum","ANSNELALAA"},
     { "Verrucosispora maris","AKHNRADFALAA"},
     { "Vesicomyosocius okutanii","ENENNYALAA"},
     { "Vibrio anguillarum","ANDENYALAA"},
     { "Vibrio campbellii","ANDENYALAA"},
     { "Vibrio cholerae","ANDENYALAA"},
     { "Vibrio Ex25","ANDENYALAA"},
     { "Vibrio fischeri","ANDENYALAA"},
     { "Vibrio furnissii","ANDENYALAA"},
     { "Vibrio parahaemolyticus","ANDENYALAA"},
     { "Vibrio parahemolyticus","ANDENYALAA"},
     { "Vibrio sp. EJY3","ANDENYALAA"},
     { "Vibrio sp. Ex25","ANDENYALAA"},
     { "Vibrio splendidus","ANDENYALAA"},
     { "Vibrio vulnificus","ANDENYALAA"},
     { "Waddlia chondrophila","ADLDLATAAVAA"},
     { "Weeksella virosa","GNEEYALAA"},
     { "Weissella koreensis","AKNSNNLAFAA"},
     { "Wigglesworthia brevipalpis","AKHKYNEPALLAA"},
     { "Wigglesworthia glossinidia","AKHKYNEPALLAA"},
     { "Wolbachi.sp","ANDNFAAEDNVDAIAA"},
     { "Wolbachia endosymbiont","ANDNFAAEEYRVAA"},
     { "Wolbachia sp. 2 (Brugi)","ANDNFAAEGDVAVAA"},
     { "Wolbachia sp. 3 (Culex)","ANDNFAAEDNVALAA"},
     { "Wolbachia sp. 4 (Dros.)","ANDNFAAEEYRVAA"},
     { "Wolinella succinogenes","ALSSHPKRGKRLGLPITSALGA"},
     { "Xanthobacter autotrophicus","ANDNYAPVAQAA"},
     { "Xanthomonas albilineans","ANDDNYALAA"},
     { "Xanthomonas axonopodis","ANDDNYGSDFAIAA"},
     { "Xanthomonas campestris 1","ANDDNYGSDFAIAA"},
     { "Xanthomonas campestris 2","ANDDNYGSDSAIAA"},
     { "Xanthomonas oryzae","ANDDNYGSDFAIAA"},
     { "Xenorhabdus bovienii","ANDENYALAA"},
     { "Xenorhabdus nematophila","ANDENYALAA"},
     { "Xylanimonas cellulosilytica","ADNTRNDFALAA"},
     { "Xylella fastidiosa 1","ANEDNFAVAA"},
     { "Xylella fastidiosa 2","ANEDNFALAA"},
     { "Xylella fastidiosa 3","ANEDNFAIAA"},
     { "Xylella fastidiosa 4","ANEDNFALAA"},
     { "Yersinia bercovieri","ANDSQYESAALAA"},
     { "Yersinia enterocolitica","ANDSQYESAALAA"},
     { "Yersinia frederiksenii","ANDENYALAA"},
     { "Yersinia intermedia","ANDSQYESAALAA"},
     { "Yersinia mollaretii","ANDSQYESAALAA"},
     { "Yersinia pestis","ANDENYALAA"},
     { "Yersinia pseudotuberculosis","ANDENYALAA"},
     { "Zobellia galactanivorans","GENNYALAA"},
     { "Zunongwangia profunda","GENNYALAA"} };

  

/* TOOLS */

char upcasec(char c)
{ return((c >= 'a')?c-32:c); }

int length(char *s)
{ int i = 0;
  while (*s++) i++;
  return(i); }

char *softmatch(char *s, char *key)
{ while (upcasec(*key) == upcasec(*s))
   { if (!*key++) return(s);
     s++; }
  if (*key) return(NULL);
  return(s); }

char *strpos(char *s, char *k)
{ char c,d;
  int i;
  d = *k;
  while (c = *s)
   { if (c == d)
       { i = 0;
         do if (!k[++i]) return(s); while (s[i] == k[i]); }
     s++; }
  return(NULL); }


char *softstrpos(char *s, char *k)
{ char c,d;
  int i;
  d = upcasec(*k);
  while (c = *s)
   { if (upcasec(c) == d)
       { i = 0;
         do if (!k[++i]) return(s);
         while (upcasec(s[i]) == upcasec(k[i])); }
     s++; }
  return(NULL); }

char *wildstrpos(char *s, char *k)
{ char c,d;
  int i;
  d = upcasec(*k);
  while (c = *s)
   { if ((upcasec(c) == d) || (d == '*'))
       { i = 0;
         do if (!k[++i]) return(s);
         while ((upcasec(s[i]) == upcasec(k[i])) || (k[i] == '*')); }
     s++; }
  return(NULL); }

char *marginstring(char *s, char *k, int margin)
{ char c,d;
  int i,j;
  j = 0;
  d = *k;
  while (c = *s)
   { if (c == d)
       { i = 0;
         do if (!k[++i]) return(s); while (s[i] == k[i]); }
     s++;
     if (++j >= margin) break; }
  return(NULL); }


int margindetect(char *line, int margin)
{ int i;
  char c,*s;
  i = 0;
  s = line;
  while (c = *s++)
   { if (!space(c)) break;
     if (c == '\t') i += 7;
     if (++i >= margin) return(0); } 
  if (c == '\n') return(0);
  if (c == '\r') return(0);
  if (c == '\0') return(0);
  return(1); }


char *backword(char *line, char *s, int n)
{
int spzone;
if (space(*s))
 { spzone = 1; }
else
 { spzone = 0;
   n++; }
while (s > line)
 { if (space(*s))
    { if (spzone == 0)
       { spzone = 1;
         if (--n <= 0) 
          return(++s); }}
   else spzone = 0;
   s--; }
if (!space(*s))
 if (n <= 1) return(s);
return(NULL);
}
     
char *dconvert(char *s, double *r)
{ static char zero='0',nine='9';
  int shift,expshift,sgn,expsgn,exponent;
  char c,limit;
  double result;
  shift = 0;
  expshift = 0;
  sgn = 1;
  expsgn = 1;
  limit = 0;
  exponent = 0;
  result = 0.0;
  if ((c = *s) == '-')
   { sgn = -1;
     c = *++s; }
  else if (c == '+') c= *++s;
  if (c >= zero)
   if (c <= nine)
    { result = (double)(c - zero);
      while ((c = *++s) >= zero)
       { if (c > nine) break;
         if (++limit < 15) result = result*10.0 + (double)(c - zero); }}
  if (c == '.')
   while ((c = *++s) >= zero)
    { if (c > nine) break;
      if (++limit < 15)
       { result = result*10.0 + (double)(c - zero);
         shift++; }}
  if ((c == 'E')||(c == 'e')||(c == 'D')||(c == 'd'))
    { if ((c = *++s) == '-')
       { expsgn = -1;
         c = *++s; }
      else
       if (c == '+') c = *++s;
      if (c >= zero)
       if (c <= nine)
         { exponent = c - zero;
           while ((c = *++s) >= zero)
            { if (c > nine) break;
              exponent = exponent*10 + c - zero;
              if (++expshift > 3) break; }}}
  result *= (double)sgn;
  exponent = exponent*expsgn - shift;
  if (exponent >= 0)
    while (exponent--) result *= 10.0;
  else
    while (exponent++) result /= 10.0;
  (*r) *= 0.01*result;
  return(s); }

char *lconvert(char *s, long *r)
{ static char zero='0',nine='9';
  long sgn;
  long result;
  char c;
  sgn = 1L;
  result = 0L;
  if ((c = *s) == '-')
   { sgn = -1L;
     c = *++s; }
  else if (c == '+') c= *++s;
  if (c >= zero)
   if (c <= nine)
    { result = (long)(c - zero);
      while ((c = *++s) >= zero)
       { if (c > nine) break;
         result = result*10L + (long)(c - zero); }}
  *r = result * sgn;
  return(s); }


char *getlong(char *line, long *l)
{ static char zero='0',nine='9';
  char c1,c2,*s;
  if (!line) return(NULL);
  s = line;
  while (c1 = *s) 
   { if (c1 >= zero)
      { if (c1 <= nine) return(lconvert(s,l)); }
     else 
      if ((c1 == '-') || (c1 == '+'))
       { c2 = s[1];
         if (c2 >= zero)
          if (c2 <= nine)
           return(lconvert(s,l)); }
     s++; }
  return(NULL); }
  
   

char *copy(char *from, char *to)
{ while (*to++ = *from++);
  return(--to);  }


char *copy2sp(char *from1, char *from2, char *to, int n)
{ 
char *s;
s = to;
while (from1 < from2)
 { *s++ = *from1++;
   if (--n <= 0)
    { do if (--s <= to) break;
      while (!space(*s));
      break; }}
*s = '\0'; 
return(s);
}


char *copy3cr(char *from, char *to, int n)
{ while (*to = *from++)
   { if (*to == DLIM)
      { *to = '\0';
        break; }
     if (--n <= 0)
      { *++to = '\0';
        break; }
     to++; }
  return(to); }

char *quotestring(char *line, char *a, int n)
{ char ch;
  while (ch = *line++) 
   if (ch == '"') 
    { while (ch = *line++) 
       { if (ch == '"') break;
         if (ch == ';') break;
         if (ch == '\n') break;
         if (ch == '\r') break;
         *a++ = ch;
         if (--n <= 0) break; }
      break; }
  *a = '\0';
  return(a); }
      

/* LIBRARY */


int fseekd(data_set *d, long fpos, long foffset)
{
if (d->bugmode)
 { fpos += foffset;
   if (fpos < 0L) fpos = 0L;
   if (fseek(d->f,0L,SEEK_SET)) return(EOF);
   d->filepointer = -1L;
   while (++d->filepointer < fpos)
    if (getc(d->f) == EOF) return(EOF);
   return(0); }
if (fseek(d->f,fpos,SEEK_SET)) return(EOF);
d->filepointer = fpos;
if (foffset != 0L)
 { if ((fpos + foffset) < 0L) foffset = -fpos;
   if (fseek(d->f,foffset,SEEK_CUR)) return(EOF);
   d->filepointer += foffset; }
return(0);
}


long ftelld(data_set *d)
{
if (d->bugmode) return(d->filepointer);
else return(ftell(d->f));
}


char fgetcd(data_set *d)
{
int ic;
if ((ic = getc(d->f)) == EOF) return(NOCHAR);
d->filepointer++;
return((char)ic);
}


char *fgetsd(data_set *d, char line[], int len)
{
int i,ic;
i = 0;
while (i < len)
 { if ((ic = getc(d->f)) == EOF) break;
   d->filepointer++;
   if (ic == '\r') continue;
   if (ic == '\n')
    { line[i++] = DLIM;
      break; }
   line[i++] = (char)ic; }
if (i < 1) return(NULL);
line[i] = '\0';
return(line);
}

int agene_position_check(data_set *d, int nagene, annotated_gene *agene)
{
int a;
long l,swap;
if ((agene->stop - agene->start) > MAXAGENELEN) 
 { swap = agene->stop;
   agene->stop = agene->start;
   agene->start = swap;
   agene->stop += d->aseqlen; }
if (agene->start > agene->stop) agene->stop += d->aseqlen;
l = agene->stop - agene->start;
if ((l < 1) || (l > MAXAGENELEN)) return(0);
if (agene->stop == d->aseqlen)
 { for (a = 0; a < nagene; a++)
    if (d->gene[a].start == agene->start)
     if (d->gene[a].genetype == agene->genetype)
      if (softmatch(d->gene[a].species,agene->species))
       return(0); }
return(1);
}



long process_sequence_heading(data_set *d, csw *sw)
{ int i,ic,nagene;
  long l,realstart;
  char line[STRLEN],c,*s,*sq,*sd;
  annotated_gene *agene,tmpagene;
  d->datatype = FASTA;
  fseekd(d,d->seqstart,d->seqstartoff);
  HEADING:           
  do if ((c = fgetcd(d)) == NOCHAR) return(-1L);
  while (space(c));
  if (c == '#')
   { if (!fgetsd(d,line,STRLENM1)) return(-1L);
     goto HEADING; }
  if (!fgetsd(d,d->seqname,STRLENM1)) return(-1L);
  if (c != '>')
   { s = d->seqname;
     if (upcasec(c) != 'L')
      { do if (!(c = *s++)) goto FNSN;
        while (upcasec(c) != 'L'); } 
     if (!(s = softmatch(s,"OCUS"))) goto FNSN;
     if (sd = softstrpos(d->seqname,"BP"))
      { sd = backword(d->seqname,sd,1);
        if (sd = getlong(sd,&l)) d->aseqlen = l; }
     s += 4;
     while (space(*s)) s++;
     sq = d->seqname;
     while (!space(*s)) *sq++ = *s++;
     d->aseqlen = 0L; 
     if (!fgetsd(d,line,STRLENM1)) return(-2L);
     if (sd = softstrpos(line,"DEFINITION"))
      { sd += 10;
        while (space(*sd)) sd++;
        *sq++ = ' ';
        copy(sd,sq);
        if (!fgetsd(d,line,STRLENM1)) return(-2L); }
     else copy(s,sq);
     for (i = 0; i < NS; i++) d->nagene[i] = 0;
     nagene = 0;
     while (!marginstring(line,"ORIGIN",10))
      { if (nagene >= NGFT) goto GBNL;
        agene = &(d->gene[nagene]);
        agene->comp = 0;
        agene->start = -1L;
        agene->stop = -1L;
        agene->antistart = -1L;
        agene->antistop = -1L;
        agene->permuted = 0;
        agene->pseudogene = 0;
        if (!(s = marginstring(line,"tRNA",10))) goto TMRNASEQ;
        agene->genetype = tRNA;
	    if (softstrpos(s,"complement")) agene->comp = 1;
        if (s = getlong(s,&l)) agene->start = l;
        if (s = getlong(s,&l)) agene->stop = l;
        copy("tRNA-???",agene->species);
        if (!fgetsd(d,line,STRLENM1)) return(-2L);
        while (!margindetect(line,10))
         { if (s = softstrpos(line,"product="))
            if (s = softstrpos(s,"tRNA-"))
             { s += 5;
               while (space(*s)) s++;
               copy3cr(s,agene->species+5,3); }
           if (s = softstrpos(line,"anticodon="))
            { s += 10;
              if (!(s = getlong(s,&l))) l = -1L;
              agene->antistart = l;
              if (!(s = getlong(s,&l))) l = -1L;
              agene->antistop = l; }
           if (softstrpos(line,"/pseudo")) agene->pseudogene = 1;
           if (!fgetsd(d,line,STRLENM1)) return(-2L); }
        if (agene_position_check(d,nagene,agene))
         { d->nagene[tRNA]++;
           nagene++; }
	    continue;
        TMRNASEQ:
        if (!(s = marginstring(line,"tmRNA",10))) goto CDSEQ;
        agene->genetype = tmRNA;
	    if (softstrpos(s,"complement")) agene->comp = 1;
        if (s = getlong(s,&l)) agene->start = l;
        if (s = getlong(s,&l)) agene->stop = l;
        copy("tmRNA",agene->species);
        if (!agene_position_check(d,nagene,agene)) goto GBNL;
        d->nagene[tmRNA]++;
        nagene++;
        if (!fgetsd(d,line,STRLENM1)) return(-2L);
        while (!margindetect(line,10))
         { if (softstrpos(line,"acceptor")) agene->permuted = 1;
           if (softstrpos(line,"/pseudo")) agene->pseudogene = 1;
           if (!fgetsd(d,line,STRLENM1)) return(-2L); }
        if (s = marginstring(line,"tmRNA",10))
         { tmpagene.comp = 0;
           tmpagene.start = -1L;
           tmpagene.stop = -1L;
           tmpagene.antistart = -1L;
           tmpagene.antistop = -1L;
           tmpagene.permuted = 0;
           tmpagene.pseudogene = 0;
	       if (softstrpos(s,"complement")) tmpagene.comp = 1;
           if (s = getlong(s,&l)) tmpagene.start = l;
           if (s = getlong(s,&l)) tmpagene.stop = l;
           if (!fgetsd(d,line,STRLENM1)) return(-2L);
           while (!margindetect(line,10))
            { if (softstrpos(line,"coding")) tmpagene.permuted = 1;
              if (softstrpos(line,"/pseudo")) tmpagene.pseudogene = 1;
              if (s = softstrpos(line,"/tag_peptide"))
               { if (s = getlong(s,&l)) tmpagene.antistart = l;
                 if (s = getlong(s,&l)) tmpagene.antistop = l; }
              if (!fgetsd(d,line,STRLENM1)) return(-2L); }
           if (agene->permuted && tmpagene.permuted)
            { agene->stop = tmpagene.stop;
              agene->antistart = tmpagene.antistart;
              agene->antistop = tmpagene.antistop;
              copy("tmRNA(Perm)",agene->species); }
           else
            { if (nagene >= NGFT) goto GBNL;
              agene = &(d->gene[nagene]);
              agene->comp = tmpagene.comp;
              agene->start = tmpagene.start;
              agene->stop = tmpagene.stop;
              agene->antistart = -1L;
              agene->antistop = -1L;
              agene->permuted = 0;
              agene->pseudogene = tmpagene.pseudogene;
              copy("tmRNA",agene->species);
              if (agene_position_check(d,nagene,agene))
               { d->nagene[tmRNA]++;
                 nagene++; }}}           
	    continue;
        CDSEQ:
        if (!(s = marginstring(line,"CDS",10)))
         if (!(s = marginstring(line,"mRNA",10))) 
          goto RRNA;
        agene->genetype = CDS;
	    if (softstrpos(s,"complement")) agene->comp = 1;
        if (s = getlong(s,&l)) agene->start = l;
        if (s = getlong(s,&l)) agene->stop = l;
        copy("???",agene->species);
        if (!fgetsd(d,line,STRLENM1)) return(-2L);
        while (!margindetect(line,10))
         { if (s = softstrpos(line,"gene="))
            { s += 5;
              quotestring(s,agene->species,SHORTSTRLENM1); }
           else if (s = softstrpos(line,"product="))
            { s += 8;
              quotestring(s,agene->species,SHORTSTRLENM1); }
           if (softstrpos(line,"/pseudo")) agene->pseudogene = 1;
           if (!fgetsd(d,line,STRLENM1)) return(-2L); }
        if (agene_position_check(d,nagene,agene))
         { d->nagene[CDS]++;
           nagene++; }
        continue;
        RRNA:
        if (!(s = marginstring(line,"rRNA",10))) goto GBNL;
        agene->genetype = rRNA;
	    if (softstrpos(s,"complement")) agene->comp = 1;
        if (s = getlong(s,&l)) agene->start = l;
        if (s = getlong(s,&l)) agene->stop = l;
        copy("???",agene->species);
        if (!fgetsd(d,line,STRLENM1)) return(-2L);
        while (!margindetect(line,10))
         { if (s = softstrpos(line,"gene="))
            { s += 5;
              quotestring(s,agene->species,SHORTSTRLENM1); }
           else if (s = softstrpos(line,"product="))
            { s += 8;
              quotestring(s,agene->species,SHORTSTRLENM1); }
           if (softstrpos(line,"/pseudo")) agene->pseudogene = 1;
           if (!fgetsd(d,line,STRLENM1)) return(-2L); }
        if (agene_position_check(d,nagene,agene))
         { d->nagene[rRNA]++;
           nagene++; }
        continue;
        GBNL:
        if (!fgetsd(d,line,STRLENM1)) return(-2L); }
     d->datatype = GENBANK;
     d->nagene[NS-1] = nagene;
     sw->annotated = 1;
     realstart = ftelld(d); }
  else
   { MH:
     realstart = ftelld(d);
     do if ((c = fgetcd(d)) == NOCHAR) return(-3L);
     while (space(c));
     if (c == '>')
      { if (!fgetsd(d,line,STRLENM1)) return(-3L);
        goto MH; }
     fseekd(d,realstart,0L); }
  s = d->seqname;
  i = 0;
  while ((c = *s) != '\0')
   { if (c == '\n') break;
     if (c == '\r') break;
     if (++i >= STRLEN) break;
     s++; }
  *s = '\0';
  return(realstart);
  FNSN:
  s = copy("Unnamed sequence ",d->seqname);
  fseekd(d,d->seqstart,d->seqstartoff);
  realstart = ftelld(d);
  if (fgetsd(d,line,STRLENM1)) copy3cr(line,s,50);
  fseekd(d,realstart,0L);
  return(realstart); }


int move_forward(data_set *d)
{ int ic;
  long nextbase;
  static int map[256] =
  { -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,NOBASE,-3,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-2,-4,-4,Adenine,AMBIG,Cytosine,AMBIG,-4,-4,Guanine,AMBIG,
    -4,-4,AMBIG,-5,AMBIG,AMBIG,-4,-4,-4,
    AMBIG,AMBIG,Thymine,Thymine,AMBIG,AMBIG,-4,
    AMBIG,-4,-4,-4,-4,INSERT,NOBASE,-4,Adenine,AMBIG,Cytosine,AMBIG,
    -4,-4,Guanine,AMBIG,-4,-4,AMBIG,-5,AMBIG,AMBIG,-4,-4,-4,
    AMBIG,AMBIG,Thymine,Thymine,AMBIG,AMBIG,-4,
    AMBIG,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4 };
  if (d->ps >= d->psmax)
   if (d->psmax > 0L)
    { fseekd(d,d->seqstart,d->seqstartoff);
      d->ps = 0L; }
  NL:
  if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
  SC:
  ic = map[ic];
  BS:
  if (ic >= Adenine)
   { d->ps++;
     return(ic); }
  if (ic == -2)
   { d->nextseq = ftelld(d);
     d->nextseqoff = -1L;
     return(TERM); }
  if (ic == -3)
   if (d->datatype == GENBANK)
   { if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
     if ((ic = map[ic]) != -3) goto BS;
     do if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
     while (space(ic));
     d->nextseq = ftelld(d);
     d->nextseqoff = -1L;
     return(TERM); }
  if (ic == -5)
   { nextbase = ftelld(d); 
     if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
     if (upcasec(ic) == 'O')
      { if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
        if (upcasec(ic) == 'C')
         { if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
           if (upcasec(ic) == 'U')
            { if ((ic = (int)fgetcd(d)) == NOCHAR) goto FAIL;
              if (upcasec(ic) == 'S')
               { d->nextseq = nextbase;
                 d->nextseqoff = -1L;
                 return(TERM); }}}}
     fseekd(d,nextbase,0L); }
  goto NL;
  FAIL:
  d->nextseq = -1L;
  d->nextseqoff = 0L;
  if (d->psmax > 0L)
   { d->ps = d->psmax;
     return(NOBASE); }
  else return(TERM); }



char cbase(int c)
{ static char base[7] = "acgt..";
  if (c < Adenine) return('#');
  if (c > NOBASE) return((char)c);
  return(base[c]); }




int seq_init(data_set *d, csw *sw)
{ long ngc;
  int ic;
  d->filepointer = 0;
  if ((d->seqstart = process_sequence_heading(d,sw)) < 0L) 
   { if (d->seqstart == -2L)
      fprintf(stderr,"ERROR - unable to read Genbank sequence %s\n",d->seqname);
     else if (d->seqstart == -2L)
      fprintf(stderr,"ERROR - unable to read fasta sequence %s\n",d->seqname);
     return(0); }
  d->seqstartoff = 0L;
  d->ps = 0L;
  d->psmax = -1L;
  ngc = 0L;
  while ((ic = move_forward(d)) >= Adenine)
   if (ic >= Cytosine)
    if (ic <= Guanine)
     ngc++;
  if ((d->psmax = d->ps) <= 0L) return(0);
  d->gc = (double)ngc/(double)d->psmax;
  fseekd(d,d->seqstart,d->seqstartoff);
  d->ps = 0L;
  return(1); }


char cpbase(int c)
{ static char base[7] = "ACGT..";
  if (c < Adenine) return('#');
  if (c > NOBASE) return((char)c);
  return(base[c]); }



char *aa(int *anticodon, csw *sw)
{ int p1,p2,p3;
  if ((p1 = *anticodon) >= AMBIG) return(ambig_aaname);
  if ((p2 = anticodon[1]) >= AMBIG) return(ambig_aaname);
  if ((p3 = anticodon[2]) >= AMBIG) return(ambig_aaname);
  return(aaname[aamap[sw->geneticcode][(p1<<4) + (p2<<2) + p3]]); }


char *translate(int *codon, csw *sw)
{ int p1,p2,p3,aa;
  if ((p1 = *codon) >= AMBIG) return(ambig_aaname);
  if ((p2 = codon[1]) >= AMBIG) return(ambig_aaname);
  if ((p3 = codon[2]) >= AMBIG) return(ambig_aaname);
  aa = aamap[sw->geneticcode][((3-p3)<<4)+((3-p2)<<2)+(3-p1)];
  if ((aa == SeC) || (aa == Pyl)) aa = Stop;
  return(aaname[aa]); }

char ltranslate(int *codon, gene *t, csw *sw)
{ int code,p1,p2,p3;
  if (t->genetype == CDS) code = t->asst;
  else code = sw->geneticcode;
  if ((p1 = *codon) >= AMBIG) return(ambig_aaname[0]);
  if ((p2 = codon[1]) >= AMBIG) return(ambig_aaname[0]);
  if ((p3 = codon[2]) >= AMBIG) return(ambig_aaname[0]);
  return(aaletter[aamap[code][((3-p3)<<4)+((3-p2)<<2)+(3-p1)]]); }


char ptranslate(int *codon, csw *sw)
{ int p1,p2,p3;
  if ((p1 = *codon) >= AMBIG) return(ambig_aaname[0]);
  if ((p2 = codon[1]) >= AMBIG) return(ambig_aaname[0]);
  if ((p3 = codon[2]) >= AMBIG) return(ambig_aaname[0]);
  return(aapolarity[aamap[sw->geneticcode][((3-p3)<<4)+((3-p2)<<2)+(3-p1)]]); }


int seqlen(gene *t)
{
return(t->nbase + t->nintron);
}


int aseqlen(data_set *d, annotated_gene *a)
{
int alen;
long astart,astop;
astart = a->start;
astop = a->stop;
if (astart > astop) astop += d->psmax;
alen = (int)(astop - astart) + 1;
return(alen);
}


double gc_content(gene *t)
{ int *s,*se;
  double ngc;
  static double score[6] = { 0.0,1.0,1.0,0.0,0.0,0.0 };
  ngc = 0.0;
  if ((t->nintron > 0) && (t->asst == 0))
   { s = t->eseq;
     se = s + t->intron;
     while (s < se) ngc += score[*s++];
     s = se + t->nintron;
     se = t->eseq + t->nbase + t->nintron;
     while (s < se) ngc += score[*s++]; }
  else
   { s = t->seq;
     se = s + t->nbase;
     while (s < se) ngc += score[*s++]; }
  return(ngc/(double)t->nbase); }


void write_seq(FILE *f, int *seq, int newline)
{ int i,c;
  i = 0;
  while ((c = *seq++) >= Adenine)
   { fputc(cbase(c),f);
     if (newline)
      if (++i >= 50)
       { fputc('\n',f);
         i = 0; }}
  if (i > 0) fputc('\n',f); }


int find_var_hairpin(gene *t)
{ int e,stem,vstem,loop,*sn,*sen,*pos1,*pos2,*sb,*se,*sc,*sd,*sf,*s;
  unsigned int c,cn,m;
  static unsigned int A[6] = { 0,0,0x100,0x400,0,0 };
  static unsigned int C[6] = { 0,0,0x400,0,0,0 };
  static unsigned int G[6] = { 0x100,0x400,0,0x200,0,0 };
  static unsigned int T[6] = { 0x400,0,0x200,0,0,0 };
  static unsigned int te[6] = { 0,0,0,0,0,0 };
  if (t->genetype != tRNA) return(0);
  if (t->var < 13) return(0);
  e = 0;
  sb = t->seq + t->astem1 + t->spacer1 + 2*t->dstem + t->dloop + 
       t->spacer2 + 2*t->cstem + t->cloop + t->nintron;
  sc = sb + 3;
  se = sb + t->var - 2;
  sf = se - 2;
  te[0] = A[*se];
  te[1] = C[*se];
  te[2] = G[*se];
  te[3] = T[*se];
  while (--se > sf)
   { te[0] = (te[0] >> 4) | A[*se];
     te[1] = (te[1] >> 4) | C[*se];
     te[2] = (te[2] >> 4) | G[*se];
     te[3] = (te[3] >> 4) | T[*se]; }
  while (se >= sc)
   { te[0] = ((te[0] >> 4) | A[*se]);
     te[1] = ((te[1] >> 4) | C[*se]);
     te[2] = ((te[2] >> 4) | G[*se]);
     te[3] = ((te[3] >> 4) | T[*se]);
     s = se - 5;
     sd = se - 7;
     m = te[*s];
     while (--s > sd) m = (m >> 4) + te[*s];
     while (s >= sb)
       {  m = (m >> 4) + te[*s];
          c = m & 0xf;
          if (c >= 9)
           { stem = 3;
             loop = (int)(se - s) - 3;
             sen = se;
             sn = s + 2;
             while (loop >= 6)
              { if ((cn = vbp[sen[-1]][sn[1]]) <= 0) break;
                c += cn;
                stem++;
                loop -= 2;
                sen--;
                sn++; }
             if (c > e)
              { e = c;
                pos1 = s;
                pos2 = sen;
                vstem = stem; }}
          s--; }
      se--; }
  if (e > 0)
   return((((int)(pos1 - sb)) << 10) + (((int)(pos2 - sb)) << 5) + vstem); 
  else
   return(0); }    



void write_to_library(FILE *f, gene *t, csw *sw)
{ int *s;
  static char trnatype[2][6] = { "tRNA","mtRNA" };
  s = t->seq + t->anticodon;
  fprintf(f,">%s",t->name);
  if (!softstrpos(t->name,"RNA"))
   switch (t->genetype)
    { case CDS:
       fprintf(f," CDS");
       break;
      case srpRNA:
       fprintf(f," srpRNA");
       break;
      case tmRNA:
       if (t->asst > 0) fprintf(f," Permuted");
       fprintf(f," tmRNA");
       break;
      case tRNA:
      default:
       t->varbp = find_var_hairpin(t);
       if (t->tstem == 0) fprintf(f," TV-loop");
       else if (t->dstem == 0) fprintf(f," D-loop");
       switch(t->cloop)
        { case 6:
           fprintf(f," %s-?""?""?(%c%c)",trnatype[sw->mtrna],
                     cbase(*s),cbase(*(s+1)));
           break;
          case 8:
           fprintf(f," %s-?""?""?(%c%c%c%c)",trnatype[sw->mtrna],
                     cbase(*s),cbase(s[1]),cbase(s[2]),cbase(s[3]));
           break;
          case 7:
          default:
           fprintf(f," %s-%s(%c%c%c)",trnatype[sw->mtrna],
                     aa(s,sw),cbase(*s),cbase(*(s+1)),cbase(*(s+2)));
           break; }
      break; }
  if (strpos(t->name,"bases)"))
   fprintf(f,"\n");
  else
   fprintf(f," (%d bases)\n",t->nbase);
  fprintf(f,"sequence =\n");
  write_seq(f,t->seq,1);
  if (*t->eseq >= Adenine)
   { fprintf(f,"extended sequence =\n");
     write_seq(f,t->eseq,1); }
  fprintf(f,"nbase = %d\n",t->nbase);
  fprintf(f,"sense = %d\n",t->comp);
  fprintf(f,"start = %ld\n",t->start);
  fprintf(f,"stop = %ld\n",t->stop);
  fprintf(f,"astem1 = %d\n",t->astem1);
  fprintf(f,"astem2 = %d\n",t->astem2);
  fprintf(f,"atail = %d\n",t->aatail);
  fprintf(f,"spacer1 = %d\n",t->spacer1);
  fprintf(f,"spacer2 = %d\n",t->spacer2);
  fprintf(f,"dstem = %d\n",t->dstem);
  fprintf(f,"dloop = %d\n",t->dloop);
  fprintf(f,"cstem = %d\n",t->cstem);
  fprintf(f,"cloop = %d\n",t->cloop);
  fprintf(f,"anticodon = %d\n",t->anticodon);
  fprintf(f,"nintron = %d\n",t->nintron);
  fprintf(f,"intron = %d\n",t->intron);
  fprintf(f,"asst = %d",t->asst);
  if (t->genetype == tmRNA)
   if (t->asst > 0) fprintf(f," permuted");
  fprintf(f,"\ntps = %d\n",t->tps);
  fprintf(f,"tpe = %d\n",t->tpe);
  fprintf(f,"var = %d\n",t->var);
  fprintf(f,"varbp = %d,%d,%d\n",((t->varbp >> 10)&0x1f),
            ((t->varbp >> 5)&0x1f),(t->varbp&0x1f));
  fprintf(f,"tstem = %d\n",t->tstem);
  fprintf(f,"tloop = %d\n",t->tloop);
  fprintf(f,"gc = %g\n\n",gc_content(t)); }



void init_tmrna(FILE *f, csw *sw)
{ int c,*s;
  s = sw->tmrna_struct;
  while ((c = *s++) != TERM) itmparam(cbase(c),f); }


      


int *make_tv(int *seq, char matrix[][MATY],
             int *x, int *y, int orient, int tv)
{ int i,px,py,stem;
  static int ux[4] = { 1,0,-1,0 };
  static int uy[4] = { 0,1,0,-1 };
  static int vx[4] = { 0,-1,0,1 };
  static int vy[4] = { 1,0,-1,0 };
  static int loopu[26][26] =
  { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
    { 1,1,1,0,0,-1,-1 },
    { 1,1,1,1,0,-1,-1,-1 },
    { 1,1,1,1,0,0,-1,-1,-1 },
    { 1,1,1,1,1,0,-1,-1,-1,-1 },
    { 1,1,1,1,1,0,0,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,0,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 } };
  static int loopv[26][26] =
  { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
    { 0,1,1,1,1,1,1 },
    { -1,1,1,1,1,1,1,1 },
    { -1,1,1,1,1,1,1,1,0 },
    { -1,0,1,1,1,1,1,1,1,0 },
    { -1,0,1,1,1,1,1,1,1,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0,0,0,0 },
    { -1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0 },
    { -1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0 },
    { -1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0 },
    { -1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 } };
  px = *x;
  py = *y;
  stem = 0;
  if (tv < 6)
   { px += ux[orient];
     py += uy[orient];
     i = 0;
     while (i < tv)
     { px += vx[orient];
       py += vy[orient];
       matrix[px][py] = cbase(*seq++);
       i++; }
     py += (6-i)*vy[orient];
     goto FN; }
  if (tv > 25)
   { if (tv % 2)
      stem = (tv - 25)/2;
     else
      stem = (tv - 24)/2;
     tv = tv - 2*stem; }
  i = 0;
  while (i < stem)
   { px += ux[orient];
     py += uy[orient];
     matrix[px][py] = cbase(*seq++);
     i++; }
  i = 0;
  while (i < tv)
  { px += ux[orient]*loopu[tv][i] + vx[orient]*loopv[tv][i];
    py += uy[orient]*loopu[tv][i] + vy[orient]*loopv[tv][i];
    matrix[px][py] = cbase(*seq++);
    i++; }
  px += ux[orient]*loopu[tv][i] + vx[orient]*loopv[tv][i];
  py += uy[orient]*loopu[tv][i] + vy[orient]*loopv[tv][i];
  i = 0;
  while (i < stem)
   { matrix[px][py] = cbase(*seq++);
     px -= (ux[orient]);
     py -= (uy[orient]);
     i++; }
  FN:
  *x = px;
  *y = py;
  return(seq); }


int base_match(char b1, char b2)
{ int i,s;
  static char base1[11] = "acgtgtagtg";
  static char base2[11] = "tgcatggatg";
  static int score[11] = { 2,2,2,2,1,1,3,3,3,3 };
  s = 0;
  for (i = 0; i < 10; i++)
   if (b1 == base1[i])
    if (b2 == base2[i])
     { s = score[i];
       break; }
  return(s); }


int *make_clover(int *seq, int b, int e, int stemlength,
                  char matrix[][MATY], int *x, int *y, int orient)
{ int i,px,py,pxb,pyb,pxe,pye,l,xlg,xlgd,ylgh,ylg;
  int *s,*se;
  static int ux[9] = { 1,0,-1,0,0,1,1,-1,-1 };
  static int uy[9] = { 0,1,0,-1,1,-1,1,1,-1 };
  static int vx[9] = { 0,-1,0,1,1,1,1,-1,-1 };
  static int vy[9] = { 1,0,-1,0,0,0,0,0,0 };
  static int loopu[18][18] =
  { { -1 }, { 0,-1 }, { 0,0,-1 }, { 0,1,-1,-1 }, { 0,1,0,-1,-1 },
    { 0,1,0,0,-1,-1 }, { 0,1,1,0,-1,-1,-1 }, { 0,1,1,0,0,-1,-1,-1 },
    { 0,1,1,1,0,-1,-1,-1,-1 }, { 0,1,1,1,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,0,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1 } };
  static int loopv[18][18] =
  { { 2 }, { 1,1 },  { 0,1,1 }, { -1,2,2,-1 }, { -1,1,1,2,-1 },
    { -1,1,1,1,1,-1 }, { -1,0,1,1,1,1,-1 }, { -1,0,1,1,1,1,0,-1 },
    { -1,0,1,1,1,1,0,0,-1 }, { -1,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 } };
  static int dloopu[18][18] =
  { { -1 }, { 0,-1 }, { 0,0,-1 }, { 0,1,-1,-1 }, { 0,1,0,-1,-1 },
    { 0,1,0,0,-1,-1 }, { 0,1,1,0,-1,-1,-1 }, { 0,1,1,0,0,-1,-1,-1 },
    { 0,1,1,0,0,0,-1,-1,-1 }, { 0,1,1,1,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,0,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1 } };
  static int dloopv[18][18] =
  { { 2 }, { 1,1 },  { 0,1,1 }, { -1,2,2,-1 }, { -1,1,1,2,-1 },
    { -1,1,1,1,1,-1 }, { -1,0,1,1,1,1,-1 }, { -1,0,1,1,1,1,0,-1 },
    { -1,0,1,1,1,1,1,-1,-1 }, { -1,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 } };
  static char bond1[5] = " +!:";
  static char bond2[5] = " +-.";
  px = *x;
  py = *y;
  s = seq + b;
  se = s + stemlength;
  while (s  < se)
  { matrix[px][py] = cbase(*s++);
    px += ux[orient];
    py += uy[orient]; }
  l = e - b - 2*stemlength;
  if (l < 0) l = 0;
  if (l < 18)
   { i = 0;
     if (orient == DOWN)
      { while (i < l)
         { px += ux[orient]*dloopu[l][i] + vx[orient]*dloopv[l][i];
           py += uy[orient]*dloopu[l][i] + vy[orient]*dloopv[l][i];
           matrix[px][py] = cbase(*s++);
           i++; }
         px += ux[orient]*dloopu[l][i] + vx[orient]*dloopv[l][i];
         py += uy[orient]*dloopu[l][i] + vy[orient]*dloopv[l][i]; }
     else
      { while (i < l)
         { px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
           py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i];
           matrix[px][py] = cbase(*s++);
           i++; }
         px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
         py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i]; }}
  else 
   { ylgh = ((l >> 2) - 2) >> 1;
     ylg = (ylgh << 1) + 2;
     xlgd = l - ylg - 2*ylgh + 1;
     xlg = (xlgd + 1) >> 1;
     pxb = px - ylgh*vx[orient];
     if ((pxb < 0) || (pxb >= MATX)) goto NOLOOP;
     pyb = py - ylgh*vy[orient];
     if ((pyb < 0) || (pyb >= MATY)) goto NOLOOP;
     pxe = px + xlg*ux[orient] + (ylg - ylgh + 1)*vx[orient];
     if ((pxe < 0) || (pxe >= MATX)) goto NOLOOP;
     pye = py + xlg*uy[orient] + (ylg - ylgh + 1)*vy[orient];
     if ((pye < 0) || (pye >= MATY)) goto NOLOOP;  
     for (i = 0; i < ylgh; i++)
      { px -= vx[orient];
        py -= vy[orient];
        matrix[px][py] = cbase(*s++); }
     for (i = 0; i < xlg; i++)
      { px += ux[orient];
        py += uy[orient];
        matrix[px][py] = cbase(*s++); }
     for (i = 1; i < ylg; i++)
      { px += vx[orient];
        py += vy[orient];
        matrix[px][py] = cbase(*s++); }
     px += vx[orient];
     py += vy[orient];
     if (!(xlgd & 1)) matrix[px][py] = cbase(*s++);
     for (i = 0; i < xlg; i++)
      { px -= ux[orient];
        py -= uy[orient];
        matrix[px][py] = cbase(*s++); }
     for (i = 1; i < ylgh; i++)
      { px -= vx[orient];
        py -= vy[orient];
        matrix[px][py] = cbase(*s++); }
     px -= (ux[orient] + vx[orient]);
     py -= (uy[orient] + vy[orient]); }
  goto STEMBOND;
  NOLOOP:
  px += ux[orient]*loopu[0][0] + vx[orient]*loopv[0][0];
  py += uy[orient]*loopu[0][0] + vy[orient]*loopv[0][0];
  STEMBOND:
  se = seq + e;
  s = se - stemlength;
  while (s  < se)
  { matrix[px][py] = cbase(*s++);
    i = base_match(matrix[px][py],
                    matrix[px - 2*vx[orient]][py - 2*vy[orient]]);
    switch(orient)
     { case  RIGHT:
       case  LEFT:  matrix[px - vx[orient]][py - vy[orient]] = bond1[i];
                    break;
       case  SLANTDR:
       case  SLANTUR:
       case  SLANTUL:
       case  SLANTDL:
       case  UPRIGHT:
       case  UP:
       case  DOWN:  matrix[px - vx[orient]][py - vy[orient]] = bond2[i];
                    break; }
    px -= ux[orient];
    py -= uy[orient]; }
  *x = px;
  *y = py;
  return(se); }





int *make_dv(int *seq, char matrix[][MATY], int dloop,
                  int orient, int *xp, int *yp)
{ int i,x,y;
  static int ux[5] = { 1,0,-1,0,0 };
  static int uy[5] = { 0,1,0,-1,1 };
  static int vx[5] = { 0,-1,0,1,1 };
  static int vy[5] = { 1,0,-1,0,0 };
  static int loopu[22][22] =
  { { -1 }, { -1,0 },
    { -1,-1,1 },
    { -1,-1,0,1 },
    { -1,-1,0,0,1 },
    { -1,-1,-1,0,1,1 },
    { -1,-1,-1,0,0,1,1 },
    { -1,-1,-1,-1,0,1,1,1 },
    { -1,-1,-1,-1,0,0,1,1,1 },
    { -1,-1,-1,-1,-1,0,1,1,1,1 },
    { -1,-1,-1,-1,-1,0,0,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,0,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,0,-1,0,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1,1,1,1 } };
  static int loopv[22][22] =
  { { -6 }, { -3,-3 },
    { -2,-2,-2 },
    { -2,-1,-1,-2 },
    { -1,-1,-2,-1,-1 },
    { -1,-1,-1,-1,-1,-1 },
    { 0,-1,-1,-1,-1,-1,-1 },
    { 0,-1,0,-1,-1,-1,-1,-1 },
    { 0,-1,0,-1,-1,-1,0,-1,-1 },
    { 0,0,-1,0,-1,-1,-1,0,-1,-1 },
    { 0,0,-1,0,-1,-1,-1,0,-1,0,-1 },
    { 0,0,0,-1,0,-1,-1,-1,0,-1,0,-1 },
    { 0,0,0,-1,0,-1,-1,-1,0,-1,0,0,-1 },
    { 0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,-1 },
    { 0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,-1 },
    { 0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,-1 },
    { 0,0,0,0,0,-1,0,-1,-1,-1,-1,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,0,-1 },
    { 0,0,0,0,0,0,-1,0,-1,-1,-1,-1,0,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,0,-1,0,-1,-1,-1,-1,0,0,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,0,0,0,-1 } };
  x = *xp;
  y = *yp;
  if ((dloop < 2) || (dloop > 21))
  { x--;
    y-= 6;
    seq += dloop;
    goto FN; }
  i = 0;
  while (i < dloop)
  { x += ux[orient]*loopu[dloop][i] + vx[orient]*loopv[dloop][i];
    y += uy[orient]*loopu[dloop][i] + vy[orient]*loopv[dloop][i];
    matrix[x][y] = cbase(*seq++);
    i++; }
  x += ux[orient]*loopu[dloop][i] + vx[orient]*loopv[dloop][i];
  y += uy[orient]*loopu[dloop][i] + vy[orient]*loopv[dloop][i];
  FN:
  *xp = x;
  *yp = y;
  return(seq); }
  
  
int *make_var(int *seq, char matrix[][MATY],
               int *x, int *y, int orient, int var, int varbp)
{ int i,b,e,p,px,py,pxf,pyf,l,stem;
  static int ux[4] = { 1,0,-1,0 };
  static int uy[4] = { 0,1,0,-1 };
  static int vx[4] = { 0,-1,0,1 };
  static int vy[4] = { 1,0,-1,0 };
  static int preu[4][5][4] =
   { { {0},{1},{1},{1},{1} },
     { {0},{1,1},{1,1},{1,1},{1,1} },
     { {0},{0,1,1},{0,1,1},{1,1,1},{1,1,1} },
     { {0},{0,0,1,1},{0,0,1,1},{0,1,1,1},{0,1,1,1} } };
  static int prev[4][5][4] =
   { { {0},{0},{0},{0},{-1} },
     { {0},{0,1},{0,0},{-1,0},{-1,0} },
     { {0},{0,1,0},{0,0,0},{-1,0,0},{-1,0,-1} },
     { {0},{0,1,0,0},{0,1,0,-1},{0,0,0,-1},{0,0,-1,-1} } };
  static int postu[4][5][4] =
   { { {0},{0},{1,-1},{0,-1,0},{1,-1,-1,0} },
     { {0},{0},{0,-1},{0,-1,0},{0,-1,-1,0} },
     { {0},{0},{0,-1},{0,-1,-1},{1,-1,-1,-1} },
     { {0},{0},{0,-1},{0,-1,-1},{1,-1,-1,-1} } };
  static int postv[4][5][4] =
   { { {0},{0},{0,1},{0,1,1},{0,1,1,1} },
     { {0},{0},{0,1},{0,1,1},{0,1,1,1} },
     { {0},{0},{0,1},{0,1,1},{0,1,1,1} },
     { {0},{0},{0,1},{0,1,1},{0,1,1,1} } };
  static int loopu[10][10] =
  { { 2 }, { 1,1 }, { 1,1,0 }, { 1,0,1,0 }, { 1,1,0,0,0 },
    { 1,1,1,-1,-1,1 }, { 1,1,1,0,0,-1,0 }, { 1,1,1,1,0,-1,-1,0 },
        { 1,1,1,1,0,0,-1,-1,0 }, { 1,1,1,1,1,0,-1,-1,-1,0 } };
  static int loopv[10][10] =
  { { 3 }, { 1,2 },  { 0,1,2 }, { 0,1,1,1 }, { -1,1,1,1,1 },
    { -1,0,1,1,1,1 }, { -1,-1,1,1,1,1,1 }, { -1,-1,0,1,1,1,1,1 },
        { -1,-1,-1,1,1,1,1,1,1 }, { -1,-1,-1,0,1,1,1,1,1,1 } };
  px = *x;
  py = *y;
  if (var < 0) var = 0;
  if (var > 30) var = 30;
  if (varbp > 0)
   { b = (varbp >> 10) & 0x1f;
     if (b > 3) goto NBP;
     stem = varbp & 0x1f;
     e = stem + ((varbp >> 5) & 0x1f);
     p = var - e;
     if (p < 1) goto NBP;
     if (p > 4) goto NBP;
     pxf = px + 2*ux[orient] + 3*vx[orient];
     pyf = py + 2*uy[orient] + 3*vy[orient];
     i = 0;
     while (i < b)
      { px += ux[orient]*preu[b][p][i] + vx[orient]*prev[b][p][i];
        py += uy[orient]*preu[b][p][i] + vy[orient]*prev[b][p][i];
        matrix[px][py] = cbase(*seq++);
        i++; }
     px += ux[orient]*preu[b][p][b] + vx[orient]*prev[b][p][b];
     py += uy[orient]*preu[b][p][b] + vy[orient]*prev[b][p][b];
     seq = make_clover(seq,0,e-b,stem,matrix,&px,&py,orient+SLANT);
     i = 0;
     while (i < p)
      { px += ux[orient]*postu[b][p][i] + vx[orient]*postv[b][p][i];
        py += uy[orient]*postu[b][p][i] + vy[orient]*postv[b][p][i];
        matrix[px][py] = cbase(*seq++);
        i++; }
     *x = pxf;
     *y = pyf;
     goto FIN;  }
  NBP:   
  if (var > 9)
   { if (var % 2) stem = (var - 7)/2;
     else stem = (var - 6)/2; }
  else stem = 0;
  l = var - 2*stem;
  i = 0;
  while (i < stem)
   { px += ux[orient] - vx[orient];
     py += uy[orient] - vy[orient];
     matrix[px][py] = cbase(*seq++);
     i++; }
  i = 0;
  while (i < l)
   { px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
     py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i];
     matrix[px][py] = cbase(*seq++);
     i++; }
  px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
  py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i];
  i = 0;
  while (i < stem)
   { matrix[px][py] = cbase(*seq++);
     px -= (ux[orient] - vx[orient]);
     py -= (uy[orient] - vy[orient]);
     i++; }
  *x = px;
  *y = py;
  FIN:
  return(seq); }





void remove_inserts(int *s1, int *s2)
{ int flag,c;
  flag = 0;
  while ((c = *s1++) != TERM)
   { if (c == INSERT)
      { flag = 1 - flag;
        continue; }
     if (flag) continue;
     *s2++ = c; }
  *s2 = TERM; }



void build_trna(gene *t, char matrix[][MATY], int x, int y, csw *sw)
{ int i,j,e,c,*seq;
  int rseq[150];
  static char bond2[5] = " +-.";
  t->varbp = find_var_hairpin(t);
  remove_inserts(t->seq,rseq);
  seq = rseq;
  i = 0;
  while (i < t->astem1)
  { matrix[x][y] = cbase(*seq++);
    y--;
    i++; }
  if (t->spacer1 > 0)
   { x--;
     if (t->spacer1 >= 3) matrix[x][y+1] = cbase(*seq++);
     matrix[x][y] = cbase(*seq++);
     y--;
     x--;
     if (t->spacer1 >= 2) matrix[x][y] = cbase(*seq++);
     if ((t->spacer2 < 2) || (t->spacer1 > 1))
      { x--;
        y--; }}
  if (t->dstem > 0)
   { e = 2*t->dstem + t->dloop;
     seq = make_clover(seq,0,e,t->dstem,matrix,&x,&y,LEFT);
     if (t->spacer2 > 1) x--;
     y--;
     if (t->spacer2 > 0) matrix[x][y] = cbase(*seq++);
     y--;
     if (t->spacer2 > 1)
      { if (t->spacer1 > 1) x++;
        matrix[x][y] = cbase(*seq++);
        if (t->spacer1 < 2) y--; }
     x++; }
  else
   seq = make_dv(seq,matrix,t->dloop,RIGHT,&x,&y);
  e = 2*t->cstem + t->cloop;
  seq = make_clover(seq,0,e,t->cstem,matrix,&x,&y,DOWN);
  if (t->tstem > 0)
   { seq = make_var(seq,matrix,&x,&y,RIGHT,t->var,t->varbp);
     e = 2*t->tstem + t->tloop;
     seq = make_clover(seq,0,e,t->tstem,matrix,&x,&y,RIGHT);
     y++; }
  else
   seq = make_tv(seq,matrix,&x,&y,RIGHT,t->tloop);
  e = t->astem2;
  i = 0;
  while (i < e)
  { if ((c = *seq++) < Adenine) break;
    matrix[x][y] = cbase(c);
    j = base_match(matrix[x][y],matrix[x - 2][y]);
    matrix[x - 1][y] = bond2[j];
    y++;
    i++; }
  i = 0;
  e = (sw->aataildisp)?ASTEM2_EXTD:t->aatail;
  j = (e < 2)?e:2;
  while (i < j)
  { if ((c = *seq++) < Adenine) break;
    matrix[x][y] = cbase(c);
    x++;
    y++;
    i++; }
  e -= j;
  i = 0;
  while (i < e)
  { if ((c = *seq++) < Adenine) break;
    matrix[x][y] = cbase(c);
    x++;
    i++; }}





void build_tmrna(gene *t, char matrix[][MATY], int x, int y, csw *sw)
{ int i,j,e,c,tarm,*seq;
  int rseq[2*MAXTMRNALEN+1];
  static char bond2[5] = " +-.";
  remove_inserts(t->eseq,rseq);
  seq = rseq + t->asst;
  i = 0;
  while (i < t->astem1)
  { matrix[x][y] = cbase(*seq++);
    y--;
    i++; }
  seq = make_dv(seq,matrix,t->dloop,RIGHT,&x,&y);
  tarm = 2*t->tstem + t->tloop;
  e = (t->asst > 0)?
      (t->cstem - t->dloop - t->astem1 - t->asst + 54):
      (2*t->cstem + t->cloop + t->nintron);
  seq = make_clover(seq,0,e,t->cstem,matrix,&x,&y,DOWN);
  seq = make_var(seq,matrix,&x,&y,RIGHT,t->var,t->varbp);
  seq = make_clover(seq,0,tarm,t->tstem,matrix,&x,&y,RIGHT);
  y++;
  e = t->astem2;
  i = 0;
  while (i < e)
  { if ((c = *seq++) == TERM) break;
    matrix[x][y] = cbase(c);
    j = base_match(matrix[x][y],matrix[x - 2][y]);
    matrix[x - 1][y] = bond2[j];
    y++;
    i++; }
  e = (sw->aataildisp)?ASTEM2_EXTD:t->aatail;
  j = (e < 2)?e:2;
  i = 0;
  while (i < j)
  { if ((c = *seq++) == TERM) break;
    matrix[x][y] = cbase(c);
    x++;
    y++;
    i++; }
  e -= j;
  i = 0;
  while (i < e)
  { if ((c = *seq++) == TERM) break;
    matrix[x][y] = cbase(c);
    x++;
    i++; } }


void init_matrix(char matrix[][MATY])
{ int i,j;
  for (i =0; i < MATY; i++)
   for (j = 0; j < MATX; j++) matrix[j][i] = ' '; }


void disp_matrix(FILE *f, char matrix[][MATY], int ylines)
{ int i,j,k;
  i = ylines;
  while (--i >= 0)
   { k = MATX;
     while (--k > 0) if (matrix[k][i] != ' ') break;
     for (j = 0; j <= k; j++) fputc(matrix[j][i],f);
     fputc('\n',f); }
  fputc('\n',f); }


void xcopy(char m[][MATY], int x, int y, char *s, int l)
{ int i;
  char c;
  i = 0;
  while (i < l)
   { if (x >= MATX) break;
     if (!(c = *s++)) break;
     m[x++][y] = c;
     i++; }}




int identify_tag(char tag[], int len, char (*thit)[50], int nt)
{ int i,n;
  char *s,*st,*sb,*sd;
  n = 0;
  st = tag + len;
  while (*--st == '*');
  for (i = 0; i < NTAG; i++)
   { s = st;
     sb = tagdatabase[i].tag;
     sd = sb;
     while (*++sd);
     while (*s-- == *--sd)
      { if (s < tag)
         { if (sd > sb) goto PAR;
           if (n >= nt) goto MANY;
           copy(tagdatabase[i].name,thit[n]);
           n++;
           break; }
        if (sd > sb) continue;
        PAR:
        if (n >= nt) goto MANY;
        s = copy(tagdatabase[i].name,thit[n]);
        copy(" (partial match)",s);
        n++;
        break; }}
  return(n);
  MANY:
  return(-1);  }



int peptide_tag(char tag[], int maxlen, gene *t, csw *sw)
{
int i,lx,*se;
se = t->eseq + t->tps;
lx = (t->tpe - t->tps + 1);
if (ltranslate(se+lx,t,sw) == '*')
 { lx += 3;
   if (ltranslate(se+lx,t,sw) == '*') lx += 3; }
lx /= 3;
if (lx > maxlen) lx = maxlen;
for (i = 0; i < lx; i++)
 { tag[i] = ltranslate(se,t,sw);
   se += 3; }
tag[i] = '\0';
return(lx);
}


void update_tmrna_tag_database(gene ts[], int nt, csw *sw)
{
int nn,i,k,c,lx;
char *sp,*se,*s;
char species[STRLEN],tag[100];
gene *t;
if (sw->tagend >= NTAGMAX) return;
for (i = 0; i < nt; i++)
 { t = ts + i;
   if (t->genetype != tmRNA) continue;
   s = t->name;
   se = NULL;
   while (*s)
    { if (*s == '|') se = s;
      s++; }
   if (!*se) continue;
   while (++se) if (space(*se)) break;
   if (!*se) continue;
   while (++se) if (!space(*se)) break;
   if (!*se) continue;
   if (softstrpos(se," sp. "))
    { if (!(sp = softstrpos(se,"two-piece")))
       if (!(sp = softstrpos(se,"tmRNA")))
        continue;
      while (space(sp[-1])) sp--;
      copy2sp(se,sp,species,49); }
   else
    { s = species;
      c = 2;
      while (*se)
       { if (space(*se))
         if (--c <= 0) break;
         *s++ = *se++; }
      *s = '\0'; }
   for (k = 0; k < sw->tagend; k++)
    if (softstrpos(tagdatabase[k].name,species)) break;
   if (k < sw->tagend) continue;
   copy(species,tagdatabase[sw->tagend].name);
   s = tag;
   lx = peptide_tag(s,50,t,sw);
   s += (lx - 1);
   while (*s == '*') s--;
   *++s = '\0';
   copy(tag,tagdatabase[sw->tagend].tag);
   if (++sw->tagend >= NTAGMAX) break; }
}

int string_compare(char *s1, char *s2)
{
int r;
char c1,c2;
r = 0;
while (c1 = *s1++)
 { if (!(c2 = *s2++)) break;
   r = (int)upcasec(c1) - (int)upcasec(c2);
   if (r != 0) break; }
return(r);
}

void report_new_tmrna_tags(csw *sw)
{
int k,n,sort[NTAGMAX];
for (n = 0; n < sw->tagend; n++)
 { k = n;
   while (--k >= 0)
    { if (string_compare(tagdatabase[n].name,tagdatabase[sort[k]].name) >= 0) break;
      sort[k+1] = sort[k]; }
   sort[++k] = n; }
fprintf(sw->f,"\ntmRNA tag database update:\n");
for (k = 0; k < sw->tagend; k++)
 { n = sort[k];
   fprintf(sw->f,"     { \"%s\",\"%s\"},\n",
         tagdatabase[n].name,tagdatabase[n].tag); }
fprintf(sw->f,"\n%d tmRNA peptide tags\n",sw->tagend);
fprintf(sw->f,"%d new tmRNA peptide tags\n\n",sw->tagend - NTAG);
}


void disp_peptide_tag(FILE *f, gene *t, csw *sw)
{ int i,lx,nm,nmh,c1,c2,c3,*s,*se;
  char tag[52],thit[21][50];
  fprintf(f,"Tag peptide at [%d,%d]\nTag sequence: ",t->tps+1,t->tpe+1);
  lx = peptide_tag(tag,50,t,sw);
  se = t->eseq + t->tps;
  s = se;
  for (i = 0; i < lx; i++)
  { if (i > 0) fputc('-',f);
    if ((c1 = *s++) >= AMBIG) continue;
    if ((c2 = *s++) >= AMBIG) continue;
    if ((c3 = *s++) >= AMBIG) continue;
    fputc(cbase(c1),f);
    fputc(cbase(c2),f);
    fputc(cbase(c3),f); }
  s = se;
  fprintf(f,"\nTag peptide:  ");
  for (i = 0; i < lx; i++)
  { fprintf(f,"%s",translate(s,sw));
    s += 3;
    if (i < (lx-1)) fputc('-',f); }
  fprintf(f,"\nTag peptide:  %s",tag);  
  if (sw->energydisp)
   { s = se;
     fprintf(f,"\nTag Polarity: ");
     for (i = 0; i < lx; i++)
     { fprintf(f,"%c",ptranslate(s,sw));
       s += 3; }}
  fputc('\n',f);
  nmh = identify_tag(tag,lx,thit,21);
  if (nmh > 0)
   { if (nmh > 1)
      { fprintf(f,"Match with tmRNA tags from:\n");
        i = 0;
        for (nm = 0; nm < nmh; nm++)
         { if (++i > 3)
            { fputc('\n',f);
              i = 1; }
            if (i > 1) fprintf(f,", ");
           fprintf(f,"%s",thit[nm]); }
        fputc('\n',f); }
     else
       fprintf(f,"Match with %s tmRNA tag\n",thit[0]); }
  else 
   if (nmh == -1)
    fprintf(f,"Match with many tmRNA tags\n");
   else
    fprintf(f,"Tag not identified\n");
  fputc('\n',f);  }


void sense_switch(int *seq1, int *seq2, int lseq)
{ int i,b;
  int *sseq,*cseq;
  sseq = seq1;
  cseq = seq2 + lseq;
  while (--cseq >= seq2)
   { b = *sseq++;
     if (b >= Adenine)
      { if (b <= Thymine)
         *cseq = Thymine - b;
        else
         { if (b <= NOBASE) *cseq = b;
           else *cseq = NOBASE; }}
     else *cseq = NOBASE; }}


double nenergy(gene *t, csw *sw)
{ double eref;
  if (t->genetype != tRNA) eref = sw->eref[t->genetype];
  else
   if (sw->mtrna)
    { if (t->dstem == 0) eref = mtRNAtthresh;
      else 
       if (t->tstem == 0) eref = mtRNAdthresh;
       else eref = mtRNAdtthresh; }
   else eref = sw->eref[tRNA];
  return(100.0*t->energy/eref); }


char *position(char *s, gene *t, csw *sw)
{ long start;
  start = t->start;
  if (sw->linear) if (start <= 0) start--;
  if (t->comp)
   sprintf(s,"c[%ld,%ld]",start,t->stop);
  else
   sprintf(s,"[%ld,%ld]",start,t->stop);
  return(s); }


void location(char *s, gene *t, csw *sw, char *m)
{ char sp[80];
  sprintf(s,"%s %s",m,position(sp,t,sw)); }

void disp_location(gene *t, csw *sw, char *m)
{ char sp[80];
  fprintf(sw->f,"%s %s\n",m,position(sp,t,sw)); }


char *name(gene *t, char *si, int proc, csw *sw)
{ int s[5],*ss,*sin,*sm,*s0,*s1,*s2,*s3,nintron;
  char *sb,*st;
  static char trnatype[2][6] = { "tRNA","mtRNA" };
  switch (t->genetype)
   { case CDS:
              sprintf(si,"CDS");
              break;
     case srpRNA:
              sprintf(si,"srpRNA");
              break;
     case tmRNA:
              if (sw->dispmatch)
               { if (t->asst > 0)
                  sprintf(si,"tmRNA(Perm)  ");
                 else sprintf(si,"tmRNA        "); }
              else
               { if (t->asst > 0)
                  sprintf(si,"tmRNA (Permuted)");
                 else sprintf(si,"tmRNA"); }
              break;
     case tRNA:
              ss = (proc?t->seq:t->ps);
              sm = ss + t->anticodon - 1;
              s0 = sm + 1;
              s1 = s0 + 1;
              s2 = s1 + 1;
              s3 = s2 + 1;
              nintron = t->nintron;
              if ((proc == 0) && (nintron > 0))
               { sin = ss + t->intron;
                 if (sm >= sin) sm += nintron;
                 if (s0 >= sin) s0 += nintron;
                 if (s1 >= sin) s1 += nintron;
                 if (s2 >= sin) s2 += nintron;
                 if (s3 >= sin) s3 += nintron; }
              s[0] = *sm;
              s[1] = *s0;
              s[2] = *s1;
              s[3] = *s2;
              s[4] = *s3;
              st = trnatype[sw->mtrna];
              sb = si;
              if (t->dstem == 0)
               { sprintf(sb,"D-loop ");
                 sb += 7; }
              if (t->tstem == 0)
               { sprintf(sb,"TV-loop ");
                 sb += 8; }
              if (t->cloop == 8)
               sprintf(sb,"%s-?(%s|%s)(%c%c%c%c)",st,
			  aa(s+1,sw),aa(s+2,sw),
                          cbase(s[1]),cbase(s[2]),cbase(s[3]),
                          cbase(s[4]));
              else if (t->cloop == 6)
               sprintf(sb,"%s-?(%s|%s)(%c%c)",st,
			  aa(s,sw),aa(s+1,sw),
                          cbase(s[1]),cbase(s[2]));
              else
               sprintf(sb,"%s-%s(%c%c%c)",st,
                          aa(s+1,sw),cbase(s[1]),cbase(s[2]),cbase(s[3]));
              break;
     default: *si = '\0';
              break; }
  return(si); }


void disp_intron(FILE *f, gene *t, csw *sw)
{ int i,c,*s,*sb,*se;
  char genename[100];
  if (t->nintron <= 0) return;
  name(t,genename,1,sw);
  fprintf(f,"Intron from %s\n",genename);
  fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n");
  sb = t->eseq + t->intron;
  s = sb;
  se = sb + t->nintron;
  i = 0;
  while (s < se)
   { if ((c = *s++) < Adenine) break;
     fputc(cbase(c),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f);
  fputc('\n',f);
  fprintf(f,"Intron Length: %d\n",t->nintron);
  fprintf(f,"Intron Insertion Position(%d-%d): ",t->intron,t->intron+1);
  s = sb - 5;
  for (i = 0; i < 5; i++) fputc(cbase(*s++),f);
  fprintf(f,"-Intron-");
  s = se;
  for (i = 0; i < 5; i++) fputc(cbase(*s++),f);
  fputc('\n',f);
  fputc('\n',f); }
  

void disp_fasta_seq(FILE *f, gene *t, int ns, int n, int nsp, int c, csw *sw)
{ int i,*s,*se;
  char genename[100],genepos[100];
  if (t->nintron > 0)
   { s = t->eseq;
     se = s + t->nbase + t->nintron; }
  else
   { s = t->seq;
     se = s + t->nbase; }
  name(t,genename,1,sw);
  position(genepos,t,sw);
  if (nsp > 0)
   { if (ns > 0) fprintf(f,">%d-%d%s%s\n",ns,n,genename,genepos);
     else fprintf(f,">%s%s\n",genename,genepos); }
  else
   { if (ns > 0) fprintf(f,">%d-%d %s %s\n",ns,n,genename,genepos);
     else fprintf(f,">%s %s\n",genename,genepos); }
  i = 0;
  while (s < se)
   { if (c) fputc(cpbase(*s++),f);
     else fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f); }


void disp_seq(FILE *f, gene *t, csw *sw)
{ int i,*s,*se;
  char genename[100];
  if (sw->seqdisp >= 3)
   { if (!sw->batch) fputc('\n',f);
     if (sw->seqdisp == 3) disp_fasta_seq(f,t,0,0,0,0,sw);
     else disp_fasta_seq(f,t,0,0,0,1,sw); }
  else
   { if (!sw->batch)
      { name(t,genename,1,sw);
        fprintf(f,"\nPrimary sequence for %s\n",genename);
        fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n"); }
     if (t->nintron > 0)
      { s = t->eseq;
        se = s + t->nbase + t->nintron; }
     else
      { s = t->seq;
        se = s + t->nbase; }
     i = 0;
     while (s < se)
      { fputc(cbase(*s++),f);
        if (++i >= 50)
         { fputc('\n',f);
           i = 0; }}
     if (i > 0) fputc('\n',f); }
  if (!sw->batch) 
   { fputc('\n',f);
     fputc('\n',f); }
}


void disp_gene_SVG(gene *t, char m[][MATY], csw *sw)
{
int i,x,y,xb,xe,yb,ye,xbm,xem,ybm,yem,xdiff,ydiff;
double xpos,ypos,xsc,ysc,fontsize,yv;
char genename[100];
FILE *f = sw->f;
xe = 0;
xb = MATX;
ye = 0;
yb = MATY;
for (y = 5; y <= 30; y++)
 { for (x = 0; x < MATX; x++)
    { if ((m[x][y] > ' ') && (m[x][y] <= '~'))
       { if (x < xb) xb = x;
         if (x > xe) xe = x;
         if (y < yb) yb = y;
         if (y > ye) ye = y; }}}
xbm = xb - 5;
if (xbm < 0) xbm = 0;
xem = xe + 5;
if (xem >= MATX) xem = MATX - 1;
if ((xb - xbm) > (xem - xe)) xbm = xb - (xem - xe);
else xem = xe + (xb - xbm);
xdiff = xem - xbm;
ybm = yb - 5;
if (ybm < 5) ybm = 5;
yem = ye + 5;
if (yem > 30) yem = 30;
if ((yb - ybm) > (yem - ye)) ybm = yb - (yem - ye);
else yem = ye + (yb - ybm);
xdiff = xem - xbm;
ydiff = yem - ybm + 4;
ysc = 10.0;
xsc = 0.1*(double)((int)(10.0*ysc*((double)xdiff/(double)ydiff) + 0.5));
fontsize = 1.4;
yv = 0.01*(double)((int)(18.5*(double)ydiff*fontsize/ysc + 0.5));
name(t,genename,1,sw);
if (!sw->batch) fprintf(f,"Scalable vector graphics (SVG) image:\n");
fprintf(f,"<svg xmlns='http://www.w3.org/2000/svg' version='1.1' ");
fprintf(f,"width='%gcm' height='%gcm' viewBox='0 0 %d %d'>\n",xsc,ysc,xdiff,ydiff);
fprintf(f,"<title>%s</title>\n",genename);
fprintf(f,"<g font-family='Courier New,Courier,monospace' font-size='%g' ",fontsize);
fprintf(f,"text-anchor='middle' fill='black' stroke='none'>\n");
i = 0;
for (y = ybm; y <= yem; y++)
 { ypos = (double)(ydiff - 2 - (y - ybm));
   for (x = xbm; x <= xem; x++)
    { if ((m[x][y] > ' ') && (m[x][y] <= '~') && (m[x][y] != '!'))
       { xpos = (double)(x - xbm);
         fprintf(f,"<text x='%g' y='%g'>%c</text>",xpos,ypos,m[x][y]);
         if (++i >= 4)
          { fputc('\n',f);
            i = 0; }}}}
if (i > 0) fputc('\n',f);
fprintf(f,"</g><g fill='none' stroke='black' stroke-width='0.075'>\n");
i = 0;
for (y = ybm; y <= yem; y++)
 { ypos = (double)(ydiff - 2 - (y - ybm));
   for (x = xbm; x <= xem; x++)
    { if ((m[x][y] == '!'))
       { xpos = (double)(x - xbm);
         fprintf(f,"<line x1='%g' y1='%g' x2='%g' y2='%g'/>",xpos,ypos,xpos,ypos-yv);
         if (++i >= 2)
          { fputc('\n',f);
            i = 0; }}}}
if (i > 0) fputc('\n',f);
fprintf(f,"</g></svg>\n");
}


void disp_trna_bracket_notation(FILE *f, gene *t, csw *sw)
{
int i,j,k,varbp,stem,ab,ae,hl,*s,*se,*sl,*sb,*sr;
char genename[100];
static int bplb[2] = { '.','(' };
static int bprb[2] = { '.',')' };
if (!sw->batch)
 { name(t,genename,1,sw);
   fprintf(f,"\nSecondary structure (bracket notation) for %s\n",genename); }
if (t->nintron > 0)
 { s = t->eseq;
   se = s + t->nbase + t->nintron; }
else
 { s = t->seq;
   se = s + t->nbase; }
sl = s;
while (sl < se) fputc(cbase(*sl++),f);
fputc('\n',f);
sl = s;
sr = se - t->aatail - 1;
for (i = 0; i < t->astem1; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
for (i = 0; i < t->spacer1; i++) fputc('s',f);
sl += t->spacer1;
sb = sl + t->dstem - 1;
sr = sb + t->dstem + t->dloop;      
for (i = 0; i < t->dstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
for (i = 0; i < t->dloop; i++) fputc('d',f);
sl += t->dloop;
for (i = 0; i < t->dstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
for (i = 0; i < t->spacer2; i++) fputc('s',f);
sl += t->spacer2;
sb = sl + t->cstem - 1;
sr = sb + t->cstem + t->cloop + t->nintron;
for (i = 0; i < t->cstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
hl = t->astem1 + t->spacer1 + 2*t->dstem + t->dloop + 
     t->spacer2 + t->cstem;
if (t->nintron > 0)
 { j = t->intron - hl;
   ab = t->anticodon - hl;
   ae = ab + t->cloop - 5;
   for (i = 0; i < j; i++) 
     if (i <= ae)
      if (i >= ab) fputc('A',f);
      else fputc('c',f);
     else fputc('c',f);
   for (i = 0; i < t->nintron; i++) fputc('i',f);
   for (i = j; i < t->cloop; i++) 
   if (i <= ae)
    if (i >= ab) fputc('A',f);
    else fputc('c',f);
   else fputc('c',f); }
else   
  { j = t->cloop - 4;
    ab = t->anticodon - hl;
    ae = t->cloop - ab - j;
    for (i = 0; i < ab; i++) fputc('c',f);
    for (i = 0; i < j; i++) fputc('A',f);
    for (i = 0; i < ae; i++) fputc('c',f); }
sl += (t->cloop + t->nintron);
for (i = 0; i < t->cstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
varbp = find_var_hairpin(t);
if (varbp > 0)
 { j = (varbp >> 10);
   k = (varbp >> 5) & 0x1f;
   stem = (varbp & 0x1f);
   sr = sl + k + stem - 1;
   sl += j;
   sb = sl + stem - 1;
   for (i = 0; i < j; i++) fputc('s',f);
   for (i = 0; i < stem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
   for (i = j+stem; i < k; i++,sl++) fputc('v',f);
   for (i = 0; i < stem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
   for (i = k+stem; i < t->var; i++,sl++) fputc('s',f); }
else
 { for (i = 0; i < t->var; i++) fputc('v',f);
   sl += t->var; }
sb = sl + t->tstem - 1;
sr = sb + t->tstem + t->tloop;
for (i = 0; i < t->tstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
for (i = 0; i < t->tloop; i++) fputc('t',f);
sl += t->tloop;
for (i = 0; i < t->tstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
sb = s + t->astem1 - 1;
for (i = 0; i < t->astem2; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
fputc('\n',f);
if (!sw->batch) fputc('\n',f);
}

       
      
void disp_tmrna_trnadomain_bracket_notation(FILE *f, gene *t, csw *sw)
{
int i,j,k,varbp,stem,*s,*sa,*sb,*sc,*sd,*se,*sf,*sl,*sr;
static int bplb[2] = { '.','(' };
static int bprb[2] = { '.',')' };
if (t->nintron <= 0) return;
if (!sw->batch)
 { fprintf(f,"Secondary structure (bracket notation) for tRNA domain:\n"); }
s = t->eseq;
se = s + t->nbase + t->nintron;
if (t->asst > 0)
 { sc = s + t->asst;
   sa = sc + t->astem1 + t->dloop + t->cstem;
   sd = s + 54;
   sf = s + t->intron; }
else
 { sc = s;
   sa = s + t->astem1 + t->dloop + t->cstem;
   sd = se - t->aatail - t->astem2 - 2*t->tstem - t->tloop - t->var - t->cstem;
   sf = se; }
sl = sc;
while (sl < sa) fputc(cbase(*sl++),f);
fputc('|',f);
sr = sd;
while (sr < sf) fputc(cbase(*sr++),f);
fputc('\n',f);
sl = sc;
sr = sf - t->aatail - 1;
for (i = 0; i < t->astem1; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
for (i = 0; i < t->spacer1; i++) fputc('s',f);
sl += t->spacer1;
sb = sl + t->dstem - 1;
sr = sb + t->dstem + t->dloop;      
for (i = 0; i < t->dstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
for (i = 0; i < t->dloop; i++) fputc('d',f);
sl += t->dloop;
for (i = 0; i < t->dstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
for (i = 0; i < t->spacer2; i++) fputc('s',f);
sl += t->spacer2;
sb = sl + t->cstem - 1;
sr = sd + t->cstem - 1;
for (i = 0; i < t->cstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
fputc('|',f);
sl = sd;
for (i = 0; i < t->cstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
varbp = find_var_hairpin(t);
if (varbp > 0)
 { j = (varbp >> 10);
   k = (varbp >> 5) & 0x1f;
   stem = (varbp & 0x1f);
   sr = sl + k + stem - 1;
   sl += j;
   sb = sl + stem - 1;
   for (i = 0; i < j; i++) fputc(' ',f);
   for (i = 0; i < stem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
   for (i = j+stem; i < k; i++,sl++) fputc('v',f);
   for (i = 0; i < stem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
   for (i = k+stem; i < t->var; i++,sl++) fputc(' ',f); }
else
 { for (i = 0; i < t->var; i++) fputc('v',f);
   sl += t->var; }
sb = sl + t->tstem - 1;
sr = sb + t->tstem + t->tloop;
for (i = 0; i < t->tstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
for (i = 0; i < t->tloop; i++) fputc('t',f);
sl += t->tloop;
for (i = 0; i < t->tstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
sb = sc + t->astem1 - 1;
for (i = 0; i < t->astem2; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
for (i = 0; i < t->aatail; i++) fputc('a',f);
fputc('\n',f);
if (!sw->batch) fputc('\n',f);
}    
       
       


void disp_tmrna_seq(FILE *f, gene *t, csw *sw)
{ int i,*s,*sb,*se;
  if (t->nintron <= 0) return;
  if (*(t->name) == '\0') fprintf(f,"tmRNA sequence\n\n");
  else fprintf(f,"tmRNA sequence in %s\n\n",t->name);
  fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n");
  sb = t->eseq;
  s = sb;
  se = sb + t->intron;
  i = 0;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->tps;
  while (s < se)
   { fputc(cpbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->tpe + 1;
  while (ltranslate(se,t,sw) == '*') se += 3;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->intron + t->nintron;
  while (s < se)
   { fputc(cpbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->nbase + t->nintron;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f);
  fprintf(f,"\n5' tRNA domain at [%d,%d]\n",
          1,t->intron);
  fprintf(f,"3' tRNA domain at [%d,%d]\n",
          t->intron+t->nintron+1,t->nbase+t->nintron);
  if (sw->secstructdisp & 2) disp_tmrna_trnadomain_bracket_notation(f,t,sw);
  fprintf(f,"Resume consensus sequence at [%d,%d]: ",t->tps - 6,t->tps + 11);
  s = t->eseq + t->tps - 7;
  for (i = 0; i < 18; i++) fputc(cbase(*s++),f);
  fputc('\n',f);
  fputc('\n',f);
  disp_peptide_tag(f,t,sw); }



void disp_tmrna_perm_seq(FILE *f, gene *t, csw *sw)
{ int i,*s,*sb,*se;
  if (t->nintron <= 0) return;
  if (*(t->name) == '\0') fprintf(f,"tmRNA Sequence\n\n");
  else fprintf(f,"tmRNA sequence in %s\n\n",t->name);
  fprintf(f,"Permuted\n");
  fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n");
  sb = t->eseq;
  s = sb;
  se = sb + 54;
  i = 0;
  while (s < se)
   { fputc(cpbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->intron;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->asst;
  while (s < se)
   { fputc(cpbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->asst + t->astem1 + t->dloop + t->cstem;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->tps;
  while (s < se)
   { fputc(cpbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->tpe + 1;
  while (ltranslate(se,t,sw) == '*') se += 3;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  se = sb + t->tpe + TMPTRAILER - 54;
  while (s <= se)
   { fputc(cpbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f);
  fprintf(f,"\n5' tRNA domain at [%d,%d]\n",
          t->asst+1,t->asst+t->astem1+t->dloop+t->cstem);
  fprintf(f,"3' tRNA domain at [%d,%d]\n",
          55,t->intron);
  if (sw->secstructdisp & 2) disp_tmrna_trnadomain_bracket_notation(f,t,sw);
  fprintf(f,"Resume consensus sequence at [%d,%d]: ",t->tps - 6,t->tps + 11);
  s = t->eseq + t->tps - 7;
  for (i = 0; i < 18; i++) fputc(cbase(*s++),f);
  fputc('\n',f);
  fputc('\n',f);
  disp_peptide_tag(f,t,sw); }


void disp_cds(FILE *f, gene *t, csw *sw)
{ int i,ncodon,*s,*se;
  char c;
  ncodon = t->nbase/3;
  if (!t->tps) ncodon--;
  fprintf(f,"\n%d codons, start = %c%c%c, stop = ",ncodon,
          cbase(t->seq[0]),cbase(t->seq[1]),cbase(t->seq[2]));
  s = t->seq + 3;
  while ((i = *s++) != TERM) fputc(cbase(i),f);
  if (t->tps) fprintf(f," incomplete");
  fprintf(f,"\n1   .   10    .   20    .   30    .   40    .   50\n");
  s = t->eseq;
  se = s;
  while (*se != TERM) se++;
  if (t->tps) se -= 3;
  i = 0;
  while (s < se)
   { c = ltranslate(s,t,sw);
     fputc(c,f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }
     s += 3; }
  if (i > 0) fputc('\n',f);
  if (sw->energydisp)
   fprintf(f,"Score = %lg\n",t->energy);
  fputc('\n',f);
  fputc('\n',f); }


int pseudogene(gene *t, csw *sw)
{
if (t->energy < sw->reportpsthresh) return(1);
if (t->genetype == tRNA)
 if (t->cloop != 7)
  return(1);
return(0);
}


void disp_gene(gene *t, char m[][MATY], csw *sw)
{ double gc;
  char stat[80];
  switch(t->genetype)
   { case tmRNA:
             build_tmrna(t,m,13,27,sw);
             xcopy(m,4,3,"tmRNA (tRNA domain)",19);
             break;
     case tRNA:
             build_trna(t,m,13,27,sw);
             name(t,stat,1,sw);
             xcopy(m,4,3,stat,length(stat));
             break; }
  location(stat,t,sw,"Sequence");
  xcopy(m,4,1,stat,length(stat));
  gc = gc_content(t);
  sprintf(stat,"%d bases, %%GC = %2.1f",t->nbase,100.0*gc);
  xcopy(m,4,2,stat,length(stat));
  if (sw->reportpseudogenes)
   if (pseudogene(t,sw))
    xcopy(m,4,4,"Possible Pseudogene",19);
  if (sw->energydisp)
   { sprintf(stat,"Score = %g\n",t->energy);
     xcopy(m,4,0,stat,length(stat)); }}


void disp_batch_trna(FILE *f, gene *t, csw *sw)
{ int ls,ps,*s,anticodon;
  char pos[50],species[50];
  char m[MATX][MATY];
  static char type[2][6] = { "tRNA","mtRNA" };
  static char asterisk[2] = { ' ','*'};
  s = t->seq + t->anticodon;
  ps = sw->reportpseudogenes?(pseudogene(t,sw)?1:0):0;
  if (sw->batchfullspecies)
   { switch(t->cloop)
      { case 6:
	       sprintf(species,"%s-?(%s|%s)%c",
                   type[sw->mtrna],aa(s-1,sw),aa(s,sw),asterisk[ps]);
	       break;
        case 8:
	       sprintf(species,"%s-?(%s|%s)%c",
                   type[sw->mtrna],aa(s,sw),aa(s+1,sw),asterisk[ps]);
	       break;
        case 7:
	    default:
           sprintf(species,"%s-%s%c",type[sw->mtrna],aa(s,sw),asterisk[ps]);
	       break; }}
  else
   { switch(t->cloop)
      { case 6:
        case 8:
	     sprintf(species,"%s-???%c",type[sw->mtrna],asterisk[ps]);
	     break;
        case 7:
	    default:
         sprintf(species,"%s-%s%c",type[sw->mtrna],aa(s,sw),asterisk[ps]);
	     break; }}
  position(pos,t,sw);
  ls = length(species);
  if (ls <= 10) fprintf(f,"%-10s%28s",species,pos);
  else if (ls <= 17) fprintf(f,"%-17s%21s",species,pos);
  else fprintf(f,"%-25s%13s",species,pos);
  if (sw->energydisp)
   { fprintf(f,"\t%5.1f",t->energy); }
  anticodon = 1 + t->anticodon;
  if (t->nintron > 0)                  
   if (t->intron <= t->anticodon)      
    anticodon += t->nintron;
  fprintf(f,"\t%-4d",anticodon);
  switch(t->cloop)
   { case 6:
          fprintf(f,"\t(%c%c) ",cbase(*s),cbase(s[1]));
	  break;
     case 8:
          fprintf(f,"\t(%c%c%c%c) ",
	          cbase(*s),cbase(s[1]),cbase(s[2]),cbase(s[3]));
	  break;
     case 7:
     default:
          fprintf(f,"\t(%c%c%c)",cbase(*s),cbase(s[1]),cbase(s[2]));
	  break; }
  if (t->nintron > 0)
   fprintf(f,"i(%d,%d)",t->intron+1,t->nintron);
  fputc('\n',f);
  if (sw->secstructdisp & 2) disp_trna_bracket_notation(f,t,sw);
  if (sw->secstructdisp & 4)
   { init_matrix(m);
     build_trna(t,m,13,27,sw);
     disp_gene_SVG(t,m,sw); }
  if (sw->seqdisp) disp_seq(f,t,sw);
}


void disp_batch_tmrna(FILE *f, gene *t, csw *sw)
{ int ps,tpe,*sb,*se;
  char pos[50];
  char m[MATX][MATY];
  static char permask[2][2][3] = 
   { {"  ","p "},{"* ","p*"} };
  ps = (t->energy < 100.0)?1:0;
  position(pos,t,sw);
  fprintf(f,"tmRNA%2s%31s",permask[(t->asst == 0)?0:1][ps],pos);
  if (sw->energydisp)
   { fprintf(f,"\t%5.1f\t",t->energy); }
  tpe = t->tpe;
  sb = t->eseq + t->tps;
  se = t->eseq + tpe + 1;
  while (ltranslate(se,t,sw) == '*')
   { se += 3;
     tpe += 3; }
  fprintf(f,"\t%d,%d\t",t->tps+1,tpe+1);
  while (sb < se)
   { fputc(ltranslate(sb,t,sw),f);
     sb += 3; }
  fputc('\n',f);
  if (sw->secstructdisp & 2) disp_tmrna_trnadomain_bracket_notation(f,t,sw);
  if (sw->secstructdisp & 4)
   { init_matrix(m);
     build_tmrna(t,m,13,27,sw);
     disp_gene_SVG(t,m,sw); }
  if (sw->seqdisp) disp_seq(f,t,sw);
}


void disp_batch_srprna(FILE *f, gene *t, csw *sw)
{ int ps,tpe,*sb,*se;
  char pos[50];
  static char asterisk[2] = { ' ','*'};
  ps = (t->energy < 100.0)?1:0;
  position(pos,t,sw);
  fprintf(f,"srpRNA%c   %25s",asterisk[ps],pos);
  if (sw->energydisp)
   { fprintf(f,"\t%5.1f",t->energy); }
  fputc('\n',f);
  if (sw->seqdisp) disp_seq(f,t,sw); }

void disp_batch_cds(FILE *f, gene *t, csw *sw)
{ int ps,tpe,*sb,*se;
  char pos[50];
  static char asterisk[2] = { ' ','*'};
  ps = (t->energy < 100.0)?1:0;
  position(pos,t,sw);
  fprintf(f,"CDS%c      %25s",asterisk[ps],pos);
  if (sw->energydisp)
   { fprintf(f,"\t%5.1f",t->energy); }
  fputc('\n',f);
  if (sw->seqdisp) disp_seq(f,t,sw); }

double vloop_stability(int *sb, int var, int *varbp)
{ int e,stem,vstem,loop,*sn,*sen,*pos1,*pos2,*se,*sc,*sd,*sf,*s;
  unsigned int c,cn,m;
  static unsigned int A[6] = { 0,0,0x100,0x400,0,0 };
  static unsigned int C[6] = { 0,0,0x400,0,0,0 };
  static unsigned int G[6] = { 0x100,0x400,0,0x200,0,0 };
  static unsigned int T[6] = { 0x400,0,0x200,0,0,0 };
  static unsigned int te[6] = { 0,0,0,0,0,0 };
  e = 0;
  sc = sb + 3;   
  se = sb + var - 2; 
  sf = se - 2;
  te[0] = A[*se];
  te[1] = C[*se];
  te[2] = G[*se];
  te[3] = T[*se];
  while (--se > sf)
   { te[0] = (te[0] >> 4) | A[*se];
     te[1] = (te[1] >> 4) | C[*se];
     te[2] = (te[2] >> 4) | G[*se];
     te[3] = (te[3] >> 4) | T[*se]; }
  while (se >= sc)
   { te[0] = ((te[0] >> 4) | A[*se]);
     te[1] = ((te[1] >> 4) | C[*se]);
     te[2] = ((te[2] >> 4) | G[*se]);
     te[3] = ((te[3] >> 4) | T[*se]);
     s = se - 5;
     sd = se - 7;
     m = te[*s];
     while (--s > sd) m = (m >> 4) + te[*s];
     while (s >= sb)
       {  m = (m >> 4) + te[*s];
          c = m & 0xf;
          if (c >= 9)
           { stem = 3;
             loop = (int)(se - s) - 3;
             sen = se;
             sn = s + 2;
             while (loop >= 6)
              { if ((cn = vbp[sen[-1]][sn[1]]) <= 0) break;
                c += cn;
                stem++;
                loop -= 2;
                sen--;
                sn++; }
             if (c > e)
              { e = c;
                pos1 = s;
                pos2 = sen;
                vstem = stem; }}
          s--; }
      se--; }
  if (e > 0)
   { *varbp = (((int)(pos1-sb))<<10) + (((int)(pos2-sb))<<5) + vstem;
     return((double)(3*(vstem - 4))); }
  else
   { *varbp = 0;
     return(-12.0); }}



double find_tag_upstream_hairpin(int *se)
{ int *sb,*sd,*sf,*sh,*s;
  unsigned int c,m,mx;
  static unsigned int A[6] = { 0,0,0,0x10000,0,0 };
  static unsigned int C[6] = { 0,0,0x10000,0,0,0 };
  static unsigned int G[6] = { 0,0x10000,0,0x10000,0,0 };
  static unsigned int T[6] = { 0x10000,0,0x10000,0,0,0 };
  static unsigned int t[6] = { 0,0,0,0,0,0 };
  mx = 0;
  sf = se - 4;
  sb = se - 20;
  t[0] = A[*se];
  t[1] = C[*se];
  t[2] = G[*se];
  t[3] = T[*se];
  while (--se > sf)
   { t[0] = (t[0] >> 4) | A[*se];
     t[1] = (t[1] >> 4) | C[*se];
     t[2] = (t[2] >> 4) | G[*se];
     t[3] = (t[3] >> 4) | T[*se]; }
  sh = se - 4;
  sd = se - 30;
  while (se > sb)
   { t[0] = ((t[0] >> 4) | A[*se]);
     t[1] = ((t[1] >> 4) | C[*se]);
     t[2] = ((t[2] >> 4) | G[*se]);
     t[3] = ((t[3] >> 4) | T[*se]);
     s = sh;
     m = t[*s];
     while (--s > sd)
       {  m = (m >> 4) + t[*s];
          c = m & 0xf;
          if (c > mx) mx = c;
          if (mx == 5) goto FND; }
     sd--;
     sh--;
     se--; }
  return(0.0);
  FND:
  return(15.0); }



double find_taghairpin(int *seq)
{ int i,*s,*sb,*se,*sf;
  unsigned int c,m,mx;
  static unsigned int A[6] = { 0,0,0,1,0,0 };
  static unsigned int C[6] = { 0,0,1,0,0,0 };
  static unsigned int G[6] = { 0,1,0,1,0,0 };
  static unsigned int T[6] = { 1,0,1,0,0,0 };
  static unsigned int t[6] = { 0,0,0,0,0,0 };
  mx = 0;
  sb = seq - 20;
  se = seq - 13;
  sf = seq - 4;
  t[0] = A[*sb];
  t[1] = C[*sb];
  t[2] = G[*sb];
  t[3] = T[*sb];
  while (++sb < se)
   { t[0] = (t[0] << 4) | A[*sb];
     t[1] = (t[1] << 4) | C[*sb];
     t[2] = (t[2] << 4) | G[*sb];
     t[3] = (t[3] << 4) | T[*sb]; }
  while (sb < sf)
   { t[0] = ((t[0] << 4) | A[*sb]) & 0xffffffff;
     t[1] = ((t[1] << 4) | C[*sb]) & 0xffffffff;
     t[2] = ((t[2] << 4) | G[*sb]) & 0xffffffff;
     t[3] = ((t[3] << 4) | T[*sb]) & 0xffffffff;
     sb++;
     s = seq + 20;
     se = seq + 2;
     m = t[*s--];
     while (s > se)
      { m = (m >> 4) + t[*s--];
        c = m & 0xf;
        if (c > mx) mx = c; }
     i = 7 - (int)mx;
     while (i-- > 0)
      { m = m >> 4;
        c = m & 0xf;
        if (c > mx) mx = c; }}
  return((double)(mx << 1)); }

double stem_energy(int *s1, int *s2, int stem)
{ int *se;
  double energy;
  static double bem[6][6] =
   { { -1.072,-0.214,-1.072, ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 3.000,-1.072, 0.000, 0.000 },
     { -1.072, 3.000,-1.072, 1.286, 0.000, 0.000 },
     {  ATBOND,-1.072, 1.286,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  se = s1 + stem;
  energy = bem[*s1++][*--s2];
  while (s1  < se)
   energy += bem[*s1++][*--s2];
  return(energy); }

double astem_energy(int *s1, int *s2, int stem)
{ int *se;
  double energy;
  static double abem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  se = s1 + stem;
  energy = abem[*s1++][*--s2];
  while (s1  < se)
   energy += abem[*s1++][*--s2];
  return(energy); }


void trna_score(FILE *f, gene *t)
{ int *s,*tpos,tarm,varbp;
  double ea,eta,evls;
  static double bem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double A[6] = { 1.0,0.0,0.0,0.0,0.0,0.0 };
  static double C[6] = { 0.0,1.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,1.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,1.0,0.0,0.0 };
  if (t->genetype != tRNA) return;
  tarm = 2*t->tstem + t->tloop;
  tpos = t->seq + t->astem1 + t->spacer1 + t->dloop + 2*t->dstem + 1 +
         2*t->cstem + t->cloop + t->var;
  s = tpos + t->tstem - 1;
  eta = 6.0*(G[s[0]] + T[s[1]] + T[s[2]] + C[s[3]]) + 3.0*A[s[1]];
  s += t->tloop - 3;
  eta += 2.0*(G[*s] + A[s[1]] + T[s[3]] + C[s[4]] + C[s[5]]);
  eta += astem_energy(tpos,tpos+tarm,t->tstem);
  eta += bem[tpos[t->tstem]][tpos[t->tstem + 4]];
  eta -= 3.0*(double)(5 - t->tstem);
  if (t->tloop > 7) eta -= 3.0*(double)(t->tloop - 7);
  else eta -= 3.0*(double)(7 - t->tloop);
  s = t->seq;
  if (t->astem1 > 7) s++;
  ea = astem_energy(s,tpos+tarm+7,7);
  if (t->var > 17) evls = vloop_stability(tpos-t->var,t->var,&varbp);
  else evls = 0.0;
  fprintf(f,"\n");
  fprintf(f,"               T-arm score: %g\n",eta);
  fprintf(f,"              A-stem score: %g\n",ea);
  fprintf(f,"          V-loop stability: %g\n",evls);
  fprintf(f,"\n"); }


void tmrna_score(FILE *f, gene *t, csw *sw)
{ int r,j,te,*s,*sb,*se,*tpos,tarm;
  double e,er,et,eal,esp,ed,ec,ea,egga,etcca,egg,eta,edgg;
  double ehairpin,euhairpin;
  static int gtem[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
  static double tagend_score[4] = { 36.0, 66.0, 62.0, 72.0 };
  static int nps[126] =
   { 0,0,0,0,
     0,0,0,0,
     0,0,0,0,
     1,1,1,1,
     0,0,0,0,
     1,1,1,1,
     0,0,0,0,
     1,1,1,1,
     0,0,0,0,
     1,1,1,1,
     1,1,1,1,
     1,1,1,1,
     2,1,2,1,
     0,0,0,0,
     2,1,1,1,
     1,1,1,1,
     0,0,0,0,
     0,0,0,0,
     0,0,0,0,
     0,0,0,0,
     0,0,0,0,
     0,0,0,0 };
  static double bem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double A[6] = { 1.0,0.0,0.0,0.0,0.0,0.0 };
  static double C[6] = { 0.0,1.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,1.0,0.0,0.0,0.0 };
  static double K[6] = { 0.0,0.0,1.0,1.0,0.0,0.0 };
  static double R[6] = { 1.0,0.0,1.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,1.0,0.0,0.0 };
  static double Y[6] = { 0.0,1.0,0.0,1.0,0.0,0.0 };
  static double nA[6] = { 0,1.0,1.0,1.0,1.0,1.0 };
  static double nV[6] = { 0,0,0,1.0,1.0,1.0 };
  static double nM[6] = { 0,0,1.0,1.0,1.0,1.0 };
  if (t->genetype != tmRNA) return;
  tarm = 2*t->tstem + t->tloop;
  s = t->eseq + t->tps - 7;
  er = A[s[1]]+2.0*T[s[2]]+C[s[2]]+3.0*A[s[3]]+R[s[4]]+Y[s[6]]+
       3.0*G[s[7]]+C[s[8]];
  if (sw->tmstrict) er -= (nV[s[10]] + nV[s[11]] + nM[s[14]] + nA[s[17]]);
  er *= 4.0;
  s = t->eseq + t->tpe - 8;
  te = ((nps[(s[0]<<4) + (s[1]<<2) + s[2]] & 1) << 1)
       | (nps[(s[3]<<4) + (s[4]<<2) + s[5]] & 1);
  et = tagend_score[te];
  if (sw->tmstrict)
   { eal = 0.0;
     j = -3;
     while (j < 6)
      { te = s[j++];
        te = (te << 2) | s[j++];
        if (te == 9) eal = (double)(11 + 2*((j + 1)/3));
        j++; }
     ehairpin = find_taghairpin(s + 8);
     euhairpin = find_tag_upstream_hairpin(t->eseq + t->tps - 10); }
  else
   { eal = 15.0;
     ehairpin = 16.0;
     euhairpin = 15.0; }
  tpos = t->eseq;
  if (t->asst > 0)
   { tpos += t->cstem + t->var + 54;
     ed = 0.0; }
  else
   { tpos += t->astem1 + t->dloop + 2*t->cstem + t->nintron + t->var;
     ed = 0.001*(double)(t->tps - (long)(tpos - t->eseq)); }
  s = tpos + t->tstem - 10;
  e = K[s[0]] + G[s[1]] + A[s[2]];
  egga = K[s[1]] + G[s[2]] + A[s[3]];
  if (e > egga) egga = e;
  egga *= 6.0;
  if (egga < 18.0) egga = 0.0;
  s = tpos + tarm + 4;
  etcca = 10.0*(T[s[0]] + C[s[1]] + C[s[2]] + A[s[3]]);
  s = t->eseq + t->asst;
  egg = 7.0*(G[s[1]] + G[s[2]]);
  edgg = 0.0;
  s = t->eseq + t->asst + t->astem1;
  sb = s + 3;
  se = s + 7;
  r = gtem[*sb++];
  while (sb < se)
   { r = (r >> 4) + gtem[*sb++];
     if ((r & 3) == 2)
      { edgg = 14.0;
        break; }}
  s = tpos + t->tstem - 1;
  if (sw->tmstrict && (t->asst == 0))
   eta = 6.0*(G[s[0]] + T[s[1]] + T[s[2]] + C[s[3]]) + 3.0*A[s[1]];
  else
   eta = 6.0*(G[s[0]] + (G[s[1]] + T[s[1]]) +
         (G[s[2]] + T[s[2]]) + C[s[3]]) + 3.0*A[s[1]];
  s += t->tloop - 3;
  eta += 2.0*(G[*s] + A[s[1]] + T[s[3]] + C[s[4]] + C[s[5]]);
  eta += astem_energy(tpos,tpos+tarm,t->tstem);
  eta += bem[tpos[t->tstem]][tpos[t->tstem + 4]];
  eta -= 3.0*(double)(5 - t->tstem);
  if (t->tloop > 7) eta -= 3.0*(double)(t->tloop - 7);
  else eta -= 3.0*(double)(7 - t->tloop);
  eta *= 1.59;
  s = t->eseq + t->asst + t->astem1 + t->dloop;
  ec = stem_energy(s,tpos-t->var,t->cstem);
  s = t->eseq + t->asst;
  ea = astem_energy(s,tpos+tarm+t->astem1,t->astem1);
  esp = ((t->tpe - t->tps) < 24)?-15.0:0.0;
  e = er + et +  ed + eal + esp + egga + egg + etcca + eta + ec + ea +
      edgg + ehairpin + euhairpin;
  fprintf(f,"\n");
  fprintf(f,"     Resume sequence score: %g\n",er);
  fprintf(f,"Resume-Tarm distance score: %g\n",ed);
  fprintf(f,"         Tag peptide score: %g\n",et);
  fprintf(f,"     Tag end alanine score: %g\n",eal);
  fprintf(f,"         Short tag penalty: %g\n",esp);
  fprintf(f,"         Tag hairpin score: %g\n",ehairpin);
  fprintf(f,"Tag upstream hairpin score: %g\n",euhairpin);
  fprintf(f,"          V-loop GGA score: %g\n",egga);
  fprintf(f,"           A-stem GG score: %g\n",egg);
  fprintf(f,"         A-stem TCCA score: %g\n",etcca);
  fprintf(f,"           D-loop GG score: %g\n",edgg);
  fprintf(f,"               T-arm score: %g\n",eta);
  fprintf(f,"              C-stem score: %g\n",ec);
  fprintf(f,"              A-stem score: %g\n",ea);
  fprintf(f,"     C-stem + A-stem score: %g\n",ea + ec);
  fprintf(f,"               Total score: %g\n",e);
  fprintf(f,"          Normalised score: %g\n",nenergy(t,sw));
  fprintf(f,"\n"); }



int find_tstems(int *s, int ls, trna_loop hit[], int nh, csw *sw)
{ int i,r,c,tstem,tloop,ithresh1;
  int *s1,*s2,*se,*ss,*si,*sb,*sc,*sf,*sl,*sx,*tem;
  double ec,energy,penalty,thresh2;
  static double bem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double A[6] = { 2.0,0.0,0.0,0.0,0.0,0.0 };
  static double C[6] = { 0.0,2.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,2.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,2.0,0.0,0.0 };
  static int tem_trna[6] =
   { 0x0100, 0x0002, 0x2000, 0x0220, 0x0000, 0x0000 };
  static int tem_tmrna[6] =
   { 0x0100, 0x0002, 0x2220, 0x0220, 0x0000, 0x0000 };
  i = 0;
  tem = (sw->tmrna || (sw->threshlevel < 1.0))?tem_tmrna:tem_trna;
  ithresh1 = (int)sw->ttscanthresh;
  thresh2 = sw->ttarmthresh;
  ss = s + sw->loffset;
  si = ss + 4 - 1;
  sl = s + ls - sw->roffset + 5 + 3;
  r = tem[*si++];
  r = (r >> 4) + tem[*si++];
  r = (r >> 4) + tem[*si++];
  while (si < sl)
   { r = (r >> 4) + tem[*si++];
     if ((c = (r & 0xF)) < ithresh1) continue;
     sb = si - 7;
     sf = sb + 13;
     ec = (double)(3*c);
     for (tstem = 4; tstem <= 5; tstem++)
      { if (sb >= (sl-8)) goto NX;
        sc = sf;
        sx = si - 2;
        for (tloop = 5; tloop <= 9; tloop++)
         { if (tloop > 7)
            penalty = 3.0*(double)(tloop - tstem - 2);
           else
            penalty = 3.0*(double)(12 - tloop - tstem);
           s1 = sb;
        s2 = sc;
           se = s1 + tstem;
           energy = ec + bem[*se][se[4]] + bem[*s1++][*--s2] - penalty;
           while (s1  < se) energy += bem[*s1++][*--s2];
           energy += G[*sx] + A[sx[1]] + T[sx[3]] + C[sx[4]] + C[sx[5]];
           if (energy >= thresh2)
            { if (i >= nh)
               { fprintf(stderr,"Too many tstem hits\n");
                 goto FN; }
              hit[i].pos = sb;
              hit[i].loop = tloop;
              hit[i].stem = tstem;
              hit[i].energy = energy;
              i++; }
           sx++;
           sc++; }
        NX:
        if (--sb < ss) break;
        sf++; }}
  FN:
  return(i); }




int find_astem5(int *si, int *sl, int *astem3, int n3,
                trna_loop hit[], int nh, csw *sw)
{ int i,k;
  int *s1,*s2,*se;
  unsigned int r,tascanthresh;
  double tastemthresh,energy;
  static unsigned int tem[6] = { 0,0,0,0,0,0 };
  static unsigned int A[6] = { 0,0,0,2,0,0 };
  static unsigned int C[6] = { 0,0,2,0,0,0 };
  static unsigned int G[6] = { 0,2,0,1,0,0 };
  static unsigned int T[6] = { 2,0,1,0,0,0 };
  static double abem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  tascanthresh = (unsigned int)sw->tascanthresh;
  tastemthresh = sw->tastemthresh;
  i = 0;
  sl += n3;
  se = astem3 + n3 - 1;
  tem[0] = A[*se];
  tem[1] = C[*se];
  tem[2] = G[*se];
  tem[3] = T[*se];
  while (--se >= astem3)
   { tem[0] = (tem[0] << 4) + A[*se];
     tem[1] = (tem[1] << 4) + C[*se];
     tem[2] = (tem[2] << 4) + G[*se];
     tem[3] = (tem[3] << 4) + T[*se]; }
  r = tem[*si++];
  k = 1;
  while (++k < n3) r = (r >> 4) + tem[*si++];
  while (si < sl)
   { r = (r >> 4) + tem[*si++];
     if ((r & 15) >= tascanthresh)
      { s1 = astem3;
        s2 = si;
        se = s1 + n3;
        energy = abem[*s1++][*--s2];
        while (s1  < se)
         energy += abem[*s1++][*--s2];
        if (energy >= tastemthresh)
         { if (i >= nh)
            { fprintf(stderr,"Too many astem5 hits\n");
              goto FN; }
           hit[i].pos = si - n3;
           hit[i].energy = energy;
           i++; }}}
  FN:
  return(i); }



/*
Resume consensus sequence is: WAUARNYGCNAANNANNA
Williams, K. P., Martindale, K. A. & Bartel, D. P.  (1999)
EMBO J. 18, 5423-5433
A more general consensus sequence is NATARNYGCNRVNNMNNH
aragorn strict search uses NATARNYGCNRVNNMNNA
aragorn relaxed search uses NATARNYGC
R = A or G
Y = C or T
W = A or T
V = A or C or G
M = A or C
H = A or C or T
K = G or T
*/

int find_resume_seq(int *s, int ls, trna_loop hit[], int nh, csw *sw)
{ int e,i,j,k,a,aa[3],*si,*sb,*sf,*st,*sl;
  double al;
  unsigned int r,c,thresh;
  static int nps[105] =
   { 0,0,0,0, 0,0,0,0,
     0,0,0,0, 1,1,1,1,
     0,0,0,0, 1,1,1,1,
     0,0,0,0, 1,1,1,1,
     0,0,0,0, 1,1,1,1,
     1,1,1,1, 1,1,1,1,
     0,1,0,1, 0,0,0,0,
     0,1,1,1, 1,1,1,1,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,0 };
  static double score[4] = { 36.0, 66.0, 62.0, 72.0 };
  static unsigned int tem[6] =
   { 0x10310000, 0x01000101, 0x00010030,
     0x02000100, 0x00000000, 0x00000000 };
  static int A[6] = { 0,1,1,1,1,1 };
  static int V[6] = { 0,0,0,1,1,1 };
  static int M[6] = { 0,0,1,1,1,1 };
  thresh = (unsigned int)sw->tmrthresh;
  i = 0;
  sl = s + ls;
  r = tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  if (sw->tmstrict)
    while (s < sl)
     { r = (r >> 4) + tem[*s++];
       if ((c = (r & 0xF)) < thresh) continue;
       c -= (V[s[1]] + V[s[2]] + M[s[5]] + A[s[8]]);
       if (c < thresh) continue;
       if (i >= nh) goto FL;
       st = s - 2;
       si = st;
       sb = st + MINTAGDIST + 2;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (*si++ != Thymine)
           si++;
          else
           if (*si == Adenine)
            { if (!(*++si & 5)) goto ST1; }
           else
            if (*si == Guanine)
             { if (*++si == Adenine) goto ST1; }
            else si++;
          si++; }
       continue;
       ST1:
       if (si < sb) continue;
       al = 0.0;
       k = 0;
       j = -11;
       while (j < -2)
     { a = si[j++];
          a = (a << 2) | si[j++];
       if (a == 9) al = (double)(11 + 2*((j + 9)/3));
          a = (a << 2) | si[j++];
          aa[k++] = a; }
       hit[i].pos = st;
       hit[i].stem = (int)(si - st);
       e = (nps[aa[1]] << 1) | (nps[aa[2]]);
       hit[i].energy = (double)(c << 2) + score[e] + al +
                       find_taghairpin(si) +
                       find_tag_upstream_hairpin(st-10);
       i++; }
  else
    while (s < sl)
     { r = (r >> 4) + tem[*s++];
       if ((c = (r & 0xF)) < thresh) continue;
       if (i >= nh) goto FL;
       st = s - 2;
       si = st + MINTAGDIST;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (*si++ != Thymine)
           si++;
          else
           if (*si == Adenine)
            { if (!(*++si & 5)) goto ST2; }
           else
            if (*si == Guanine)
             { if (*++si == Adenine) goto ST2; }
            else si++;
          si++; }
       continue;
       ST2:
       hit[i].pos = st;
       hit[i].stem = (int)(si - st);
       e = (nps[(si[-8] << 4) | (si[-7] << 2) | si[-6]] << 1) |
           (nps[(si[-5] << 4) | (si[-4] << 2) | si[-3]]);
       hit[i].energy = 46.0 + (double)(c << 2) + score[e];
       i++; }
  FN:
  return(i);
  FL:
  fprintf(stderr,"Too many resume sequence hits\n");
  goto FN; }


int *base_copy3(int *from, int *to, int n)
{ while (n-- > 0) *to++ = *from++;
  *to = TERM;
  return(to);  }



void remove_intron(int *s1, int *s2, int nbase, int intron, int nintron)
{ int *s1e;
  s1e = s1 + intron;
  nbase -= intron;
  while (s1 < s1e) *s2++ = *s1++;
  s1 += nintron;
  s1e = s1 + nbase;
  while (s1 < s1e) *s2++ = *s1++;
  *s2 = TERM; }




gene *nearest_trna_gene(data_set *d, int nt, gene *t, csw *sw)
{ int n,i,comp,mtrna,mtcompov,maxintronlen,ilength;
  long a,b,c,e,score,thresh,psmax;
  static long proximity = 7*MINCTRNALEN/10;
  double energy;
  psmax = d->psmax;
  comp = t->comp;
  mtrna = sw->mtrna;
  mtcompov = sw->mtcompov;
  maxintronlen = sw->maxintronlen;
  n = -1;
  energy = INACTIVE;
  a = t->start;
  b = t->stop;
  thresh = b-a;
  if (b < a)
   { b += psmax;
     thresh += psmax;
     for (i = 0; i < nt; i++)
      { c = ts[i].start;
        e = ts[i].stop;
        if (e < c)
         { e += psmax;
           if (a > e) goto NXTW;
           if (b < c) goto NXTW;
           if (ts[i].genetype != tRNA) continue;
           if (ts[i].comp != comp) 
            { if (!mtrna) continue;
              if (mtcompov) continue; }
           if (maxintronlen > 0)
            { ilength = e - c;
              if ((2*thresh) > (5*ilength)) continue;
              if ((2*ilength) > (5*thresh)) continue; }
           score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
           if (score >= proximity)
            if (ts[i].energy < energy)
              { n = i;
                energy = ts[i].energy; }
           NXTW:
           c -= psmax;
           e -= psmax; }
        if (a > e) continue;
        if (b < c) continue;
        if (ts[i].genetype != tRNA) continue;
        if (ts[i].comp != comp)
         { if (!mtrna) continue;
           if (mtcompov) continue; }
        if (maxintronlen > 0)
         { ilength = e - c;
           if ((2*thresh) > (5*ilength)) continue;
           if ((2*ilength) > (5*thresh)) continue; }
        score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
        if (score >= proximity)
            if (ts[i].energy < energy)
              { n = i;
                energy = ts[i].energy; } }
     a -= psmax;
     b -= psmax; }
  for (i = 0; i < nt; i++)
   { c = ts[i].start;
     e = ts[i].stop;
     if (e < c)
      { e += psmax;
        if (a > e) goto NXTN;
        if (b < c) goto NXTN;
        if (ts[i].genetype != tRNA) continue;
        if (ts[i].comp != comp)
         { if (!mtrna) continue;
           if (mtcompov) continue; }
        if (maxintronlen > 0)
         { ilength = e - c;
           if ((2*thresh) > (5*ilength)) continue;
           if ((2*ilength) > (5*thresh)) continue; }
        score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
        if (score >= proximity)
            if (ts[i].energy < energy)
              { n = i;
                energy = ts[i].energy; }
        NXTN:
        c -= psmax;
        e -= psmax; }
     if (a > e) continue;
     if (b < c) continue;
     if (ts[i].genetype != tRNA) continue;
     if (ts[i].comp != comp)
      { if (!mtrna) continue;
        if (mtcompov) continue; }
     if (maxintronlen > 0)
      { ilength = e - c;
        if ((2*thresh) > (5*ilength)) continue;
        if ((2*ilength) > (5*thresh)) continue; }
     score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
     if (score >= proximity)
      if (ts[i].energy < energy)
       { n = i;
         energy = ts[i].energy; } }
  if (n >= 0) return(ts + n);
  return(NULL); }


gene *nearest_tmrna_gene(data_set *d, int nt, gene *t)
{ int n,i,comp;
  long a,b,c,e,score,smax,thresh,psmax;
  psmax = d->psmax;
  comp = t->comp;
  smax = -1;
  n = -1;
  a = t->start;
  b = t->stop;
  thresh = b-a;
  if (b < a)
   { b += psmax;
     thresh += psmax;
     for (i = 0; i < nt; i++)
      { c = ts[i].start;
        e = ts[i].stop;
        if (e < c)
         { e += psmax;
           if (a > e) goto NXTW;
           if (b < c) goto NXTW;
           if (ts[i].genetype != tmRNA) continue;
           if (ts[i].comp != comp) continue;
           score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
           if (score >= smax)
            if (score > smax)
             { n = i;
               smax = score; }
            else
             if (ts[i].energy < ts[n].energy)
               n = i;
           NXTW:
           c -= psmax;
           e -= psmax; }
        if (a > e) continue;
        if (b < c) continue;
        if (ts[i].genetype != tmRNA) continue;
        if (ts[i].comp != comp) continue;
        score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
        if (score >= smax)
         if (score > smax)
          { n = i;
            smax = score; }
         else
          if (ts[i].energy < ts[n].energy)
           n = i; }
     a -= psmax;
     b -= psmax; }
  for (i = 0; i < nt; i++)
   { c = ts[i].start;
     e = ts[i].stop;
     if (e < c)
      { e += psmax;
        if (a > e) goto NXTN;
        if (b < c) goto NXTN;
        if (ts[i].genetype != tmRNA) continue;
        if (ts[i].comp != comp) continue;
        score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
        if (score >= smax)
         if (score > smax)
          { n = i;
            smax = score; }
         else
          if (ts[i].energy < ts[n].energy)
           n = i;
        NXTN:
        c -= psmax;
        e -= psmax; }
     if (a > e) continue;
     if (b < c) continue;
     if (ts[i].genetype != tmRNA) continue;
     if (ts[i].comp != comp) continue;
     score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
     if (score >= smax)
      if (score > smax)
       { n = i;
         smax = score; }
      else
       if (ts[i].energy < ts[n].energy)
         n = i; }
  if ((10*smax) > (9*thresh)) return(ts + n);
  return(NULL); }


void overlap(data_set *d, int sort[], int n, int it, csw *sw)
{ int i,j,flag,cross,crosstoo;
  long a,b,e,f,a2,b2,e2,f2,psmax;
  char sname[100],s[100];
  flag = 0;
  cross = 0;
  psmax = d->psmax;
  a = ts[it].start;
  b = ts[it].stop;
  if (b < a)
   { a2 = a - psmax;
     b2 = b;
     b += psmax;
     cross = 1; }
  j = -1;
  while (++j < n)
   { i = sort[j];
     if (i == it) continue;
     e = ts[i].start;
     f = ts[i].stop;
     crosstoo = 0;
     if (f < e)
      { e2 = e - psmax;
        f2 = f;
        f += psmax;
        crosstoo = 1; }
     if (a <= f)
      if (b >= e)
       goto OV;
     if (crosstoo)
      if (a <= f2)
       if (b >= e2)
        goto OV;
     if (cross)
      { if (a2 <= f)
         if (b2 >= e)
          goto OV;
        if (crosstoo)
         if (a2 <= f2)
          if (b2 >= e2)
           goto OV; }
     continue;
     OV:
     if (!flag) fputc('\n',sw->f);
     name(ts+i,sname,1,sw);
     location(s,ts+i,sw,sname);
     fprintf(sw->f,"Overlap with %d: %s\n", j+1,s);
     flag = 1; }
 if (flag) fputc('\n',sw->f); }



void init_gene(int nstart, int nstop)
{ int i;
  for (i = nstart; i < nstop; i++)
   { ts[i].energy = -1.0;
     ts[i].genetype = noGENE;
     ts[i].tps = 0;
     *(ts[i].name) = '\0'; }}



gene *find_slot(data_set *d, gene *t, int *nts, csw *sw)
{ int i,newspace;
  char s1[80],s2[80],s3[80],s4[80];
  gene *tn,*tsn;
  if (sw->comp)
   { t->stop = sw->start - t->start - 1;
     t->start = t->stop - t->nbase - t->nintron + 1;
     t->comp = 1; }
  else
   { t->start += sw->start;
     t->stop = t->start + t->nbase + t->nintron - 1;
     t->comp = 0; }
  if (!sw->linear)
   { t->start = sq(t->start);
     t->stop = sq(t->stop); }
  if (t->genetype == tRNA)
   tn = nearest_trna_gene(d,*nts,t,sw);
  else
   if (t->genetype == tmRNA)  
    tn = nearest_tmrna_gene(d,*nts,t);
   else tn = NULL;
  if (tn)
   { if (t->energy <= tn->energy) return(NULL);
     copy(tn->name,t->name);
     if (sw->verbose)
      { fprintf(stderr,"%s %s ",name(t,s1,0,sw),position(s3,t,sw));
        if (sw->energydisp) fprintf(stderr,"(%g) ",nenergy(t,sw));
        fprintf(stderr,"replacing %s %s",name(tn,s2,1,sw),
        position(s4,tn,sw));
        if (sw->energydisp) fprintf(stderr," (%g)",nenergy(tn,sw));
        fprintf(stderr,"\n"); }}
  else
      { if (*nts >= sw->genespace)
         { newspace = (d->ps > 0)?(sw->genespace*(1 + d->psmax/d->ps)):
                                     (sw->genespace + NT);
           tsn = (gene *)realloc((void *)ts,newspace*sizeof(gene));
           if (tsn == NULL)
            { fprintf(stderr,"No more memory to store detected genes\n");
              fprintf(stderr,"Gene lost\n");
              return(NULL); }
           if (sw->verbose)
            fprintf(stderr,
                    "Expanding detected gene store from %d genes to %d genes\n",
                     sw->genespace,newspace);
           ts = tsn;
           init_gene(sw->genespace,newspace);
           sw->genespace = newspace; }
        copy3cr(d->seqname,t->name,99);
        tn = ts + (*nts);
        *nts = (*nts) + 1;
        if (sw->verbose)
         { fprintf(stderr,"%s at %s",name(t,s1,0,sw),position(s2,t,sw));
           if (sw->energydisp) fprintf(stderr," (%g)",nenergy(t,sw));
           fprintf(stderr,"\n"); }}
  return(tn); }


int aatail(int *s, int *ext, csw *sw)
{ int score,e;
  static int A[6] = { 1,0,0,0,0,0 };
  static int C[6] = { 0,1,0,0,0,0 };
  if (sw->aataildiv)
   { score = 0;
     e = 0;
     if (A[s[3]])
      { score++;
        e = 3; }
     if (C[s[2]])
      { score++;
        if (!e) e = 2; }
     if (C[s[1]])
      { score++;
        if (!e) e = 1; }
     if (score < 2)
      if (A[*s]) score++;
     *ext = ++e;
     return(score); }
  else
   { score = 1;
     e = 1;
     if (C[s[1]])
      { score++;
        e = 2;
        if (C[s[2]])
         { score++;
           e = 3;
           if (A[s[3]])
            { score++;
              e = 4; }}}
     *ext = e;
     return(score); }}









int find_mt_trna(data_set *d, int *seq, int lseq, int nts, csw *sw)
{ int nah,ndh,nch,nth,ncdsh,h,i,j,k,n,p,y,av,gcc,cgcc,catc,athresh;
  int igc,nbase,b8,b9,b48,b57,nc,na,nt,nti,nd,ndi,dposmap[32];
  int dl,tl,extastem,astem8,astem8d,ti,di,ser,tastem,tastem8,tastem8d;
  int astem,asteme,as,as8,aext,aext8,nbasefext,cloop,dloop,tloop,tc;
  int carm,cstem,darm,dstem,tarm,tstem,var,varbp,spacer1,spacer2,anticodon;
  int ds,dstemmotif,cloop7,mtxdetect,incds;
  int *s,*sl,*s1,*s2,*s4,*sa,*sb,*sc,*se,*sf,*sg,*si;
  int *slm,*slm1,*sle,*slb,*sld,*sge;
  int *dpos,*cpos,*cend,*tpos,*tend,*apos1,*apos2,*aend1,*aend2;
  int *clooppos,*cloopend;
  unsigned int bondtype,abondtype,mabondtype,acbondtype,cbondtype;
  unsigned int agcat,cgcat,tgcat,dbondtype,dtbondtype,tbondtype;
  unsigned int r,ct[6],cm,cv,q,tendmap[63]; 
  double gcv,e,ec,ea,eas,ed,et,ev,energy,stem_energy;
  double darmthresh,tarmthresh,tthresh,dthresh,dtthresh,thresh;
  mt_trna_cloop chit[6];
  static mt_trna_loop dhit[mtND+1];
  static mt_trna_tloop thit[mtNTH+1];
  static mt_trna_astem ahit[mtNA+1];
  static mt_cds cdshit[mtNCDS];
  gene *tn;
  static gene te =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,2,1,4,7,5,7,0,0,0,5,0,5,7,
     tRNA,0.0,0,0,0 };
  static int cAI[6] = { 8,0,0,0,8,0 };
  static int cfCI[6] = { 0,16,0,0,16,0 };
  static int cRI[6] = { 8,0,4,0,8,0 };
  static int cTI[6] = { 0,0,0,16,16,0 };
  static int cYI[6] = { 0,8,0,4,8,0 };
  static int AI[6] = { 1,0,0,0,1,0 };
  static int CI[6] = { 0,1,0,0,1,0 };
  static int GI[6] = { 0,0,1,0,1,0 };
  static int TI[6] = { 0,0,0,1,1,0 };
  static int RI[6] = { 1,0,1,0,1,0 };
  static int YI[6] = { 0,1,0,1,1,0 };
  static int WI[6] = { 1,0,0,1,1,0 };
  static unsigned int tem[6] = { 0,0,0,0,0,0 };
  static unsigned int At[6] = { 0,0,0,1,1,0 };
  static unsigned int Ct[6] = { 0,0,1,0,1,0 };
  static unsigned int Gt[6] = { 0,1,0,1,1,0 };
  static unsigned int Tt[6] = { 1,0,1,0,1,0 };
  static unsigned int cAt[6] = { 0,0,0,2,2,0 };
  static unsigned int cCt[6] = { 0,0,2,0,2,0 };
  static unsigned int cGt[6] = { 0,2,0,1,2,0 };
  static unsigned int cTt[6] = { 2,0,1,0,2,0 };
  static unsigned int aAt[6] = { 0,0,1,2,2,0 };
  static unsigned int aCt[6] = { 0,0,2,0,2,0 };
  static unsigned int aGt[6] = { 1,2,0,1,2,0 };
  static unsigned int aTt[6] = { 2,0,1,1,2,0 };
  static unsigned int dAt[6] = { 0,0,1,2,2,0 };
  static unsigned int dCt[6] = { 0,0,2,0,2,0 };
  static unsigned int dGt[6] = { 1,2,0,2,2,0 };
  static unsigned int dTt[6] = { 2,0,2,1,2,0 };
  static unsigned int clmotif[mtNCLM] =
   { 0x1321300,0x3321300,0x1323002 };
  static int dloopi[mt_DRLmaxlength+1][4] =
   { { -1 }, { -1 }, { -1 }, { -1 }, { -1 }, { -1 }, { -1 },
     { 0,2,-1 }, { 0,2,-1 }, { 0,2,3,-1 }, { 0,3,-1 }, { 0,3,-1 },
     { 0,3,4,-1 }, { 0,4,-1 }, { 0,5,-1 }, { 0,5,6,-1 }, { 0,5,6,-1 } };
  static int tloopa[12][4] =
   { { -1 }, { -1 }, { -1 }, { 0,1,-1 }, { 0,2,1,-1 }, { 4,3,2,-1 },
     { 4,3,-1 }, { 4,3,-1 }, { 4,3,-1 }, { 5,4,3,-1 }, { 5,4,-1 }, { 5,-1 } };
  static double dA[6] = { 1.0,0.0,0.0,0.0,1.0,0.0 };
  static double dT[6] = { 0.0,0.0,0.0,1.0,1.0,0.0 };
  static double C[6] = { 0.0,1.0,0.0,0.0,1.0,0.0 };
  static double G[6] = { 0.0,0.0,1.0,0.0,1.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,1.0,1.0,0.0 };
  static double AX[6] = { 0.0,-1.0,-1.0,-1.0,0.0,-1.0 };
  static double AX37[6] = { 0.0,-4.0,-1.0,-4.0,0.0,-4.0 };
  static double AXX[6] = { 0.0,-3.0,-1.5,-3.0,0.0,-3.0 };
  static double AXX37[6] = { 0.0,-4.0,-4.0,-4.0,0.0,-4.0 };
  static double AX7[6] = { 0.0,-0.7,-0.7,-0.7,0.0,-0.7 };
  static double CX[6] = { -2.0,0.0,-2.0,-1.0,0.0,-2.0 };
  static double CXX[6] = { -4.0,0.0,-4.0,-2.0,0.0,-4.0 };
  static double CX7[6] = { -0.7,0.0,-0.7,-0.7,0.0,-0.7 };
  static double TX[6] = { -1.0,-1.0,-1.0,0.0,0.0,-1.0 };
  static double TXX[6] = { -2.0,-2.0,-2.0,0.0,0.0,-2.0 };
  static double YX[6] = { -1.0,0.0,-1.0,0.0,0.0,-1.0 };
  static double tC[6] = { 0.0,0.01,0.0,0.0,0.01,0.0 };
  static double tG[6] = { 0.0,0.0,0.01,0.0,0.01,0.0 };
  static double tT[6] = { 0.0,0.0,0.0,0.01,0.01,0.0 };
  static double cA[6] = { 0.8,0.0,0.0,0.0,0.8,0.0 };
  static double cfC[6] = { 0.0,2.6,0.0,0.0,2.6,0.0 };
  static double cR[6] = { 0.8,-2.0,0.8,-0.8,0.8,-0.8 };
  static double cT[6] = { -0.8,0.0,-0.8,2.6,2.6,-0.8 };
  static double cY[6] = { -0.8,0.8,-0.8,0.8,0.8,-0.8 };
  static double loop_stab[41] =
  { 10.0,2.0,1.0,0.4,0.3,0.2,0.1,0.0,0.1,0.2,0.3,0.4,0.5,1.6,1.7,1.8,
    1.9,2.0,2.1,2.2,2.3,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
    5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7 };
  static double bem[6][6] =
   { {  mtNOBOND, mtNOBOND, mtGABOND, mtATBOND, mtATBOND, mtNOBOND },
     {  mtNOBOND, mtNOBOND, mtGCBOND, mtNOBOND, mtGCBOND, mtNOBOND },
     {  mtGABOND, mtGCBOND, mtGGBOND, mtGTBOND, mtGCBOND, mtNOBOND },
     {  mtATBOND, mtNOBOND, mtGTBOND, mtTTBOND, mtATBOND, mtNOBOND },
     {  mtATBOND, mtGCBOND, mtGCBOND, mtATBOND, mtGCBOND, mtNOBOND },
     {  mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND } };
  static double hbem[5][5] =
   { {  0.0,0.0,0.0,mtBONDSTAB+0.5*mtATBOND,mtBONDSTAB+0.5*mtATBOND },
     {  0.0,0.0,mtBONDSTAB+0.5*mtGCBOND,0.0,mtBONDSTAB+0.5*mtGCBOND },
     {  0.0,mtBONDSTAB+0.5*mtGCBOND,0.0,mtBONDSTAB+0.5*mtGTBOND,
                                            mtBONDSTAB+0.5*mtGCBOND },
     {  mtBONDSTAB+0.5*mtATBOND,0.0,mtBONDSTAB+0.5*mtGTBOND,0.0,
                                            mtBONDSTAB+0.5*mtATBOND },
     {  mtBONDSTAB+0.5*mtATBOND,mtBONDSTAB+0.5*mtGCBOND,
                                            mtBONDSTAB+0.5*mtGCBOND,
					    mtBONDSTAB+0.5*mtATBOND,
					    mtBONDSTAB+0.5*mtGCBOND } };
  tarmthresh = sw->mttarmthresh;
  tthresh = sw->mttthresh;
  dthresh = sw->mtdthresh;
  dtthresh = sw->mtdtthresh;
  ds = sw->discrim;
  extastem = sw->extastem;
  cloop7 = sw->cloop7;
  mtxdetect = sw->mtxdetect;

  /* find coding sequences */

  ncdsh = 0;

  /* find cstems */

  sc = seq + sw->loffset;
  sl = seq + lseq - sw->roffset;
  h = sc[16];
  p = sc[15];
  j = sc[14];
  k = sc[13];
  n = sc[12];
  y = sc[11];
  ct[0] = cAt[h]|(cAt[p]<<4)|(cAt[j]<<8)|(cAt[k]<<12)|
	  (cAt[n]<<16)|(cAt[y]<<20);
  ct[1] = cCt[h]|(cCt[p]<<4)|(cCt[j]<<8)|(cCt[k]<<12)|
	  (cCt[n]<<16)|(cCt[y]<<20);
  ct[2] = cGt[h]|(cGt[p]<<4)|(cGt[j]<<8)|(cGt[k]<<12)|
	  (cGt[n]<<16)|(cGt[y]<<20);
  ct[3] = cTt[h]|(cTt[p]<<4)|(cTt[j]<<8)|(cTt[k]<<12)|
	  (cTt[n]<<16)|(cTt[y]<<20);
  ct[4] = 0;
  ct[5] = 0;
  for (; sc < sl; sc++)
   { p = sc[17];
     ct[0] = (ct[0] << 4) | cAt[p];
     ct[1] = (ct[1] << 4) | cCt[p];
     ct[2] = (ct[2] << 4) | cGt[p];
     ct[3] = (ct[3] << 4) | cTt[p];
     cm = (ct[sc[4]] >> 16) + (ct[sc[3]] >> 12) + (ct[sc[2]] >> 8) +
	     (ct[sc[1]] >> 4) + ct[*sc];

   /* 7 base cloop */

     cv = (cm & 0xf0);
     athresh = 12;
     nch = 0;

  /* exclude the following cloops */ 
  /* RRnnnNN, NRnnnYN */
  /* NRnnnNN with cstem < 3 Watson-Crick basepairs or equivalent */
  /* RYnnnYN */
  /* NYnnnNN with cstem < 1 Watson-Crick basepair or equivalent */
  /* NYnnnNN with cstem < 2 Watson-Crick basepairs or equivalent */
  /* unless cloop = CTnnnAN */

     if (RI[sc[6]])
      { if (RI[sc[5]]) goto CLOOP6;
        if (YI[sc[10]]) goto CLOOP6;
        if (cv < 0x60) goto CLOOP6; }
     else
      { if (YI[sc[10]])
         if (RI[sc[5]]) goto CLOOP6;
        if (cv < 0x40)
         { if (cv < 0x20) goto CLOOP6;
           if (sc[5] != Cytosine) goto CLOOP6;
           if (sc[6] != Thymine) goto CLOOP6;
           if (sc[10] != Adenine) goto CLOOP6;
            athresh = 11; }
        else if (cv < 0x70)
         { athresh = 11;
           k = cYI[sc[5]] + cTI[sc[6]] + cRI[sc[10]] + cAI[sc[11]];
           if (sc[6] == Cytosine)
            if (sc[5] == Cytosine)
             k += 16;
            else
             if (sc[5] == Thymine)
              if (sc[11] == Adenine)
               k += 16;
           if (cv == 0x40)
            { if (k < 40) goto CLOOP6;  }
           else 
            if (cv == 0x50)
             { if (k < 28) goto CLOOP6; }
            else
             { if (k < 20) goto CLOOP6;
               athresh = 9; }}
        else
         athresh = (cv < 10)?9:8; }
     chit[0].pos = sc;
     chit[0].stem = 5;
     chit[0].loop = 7;
     chit[0].looppos = sc + 5;
     chit[0].arm = 17;
     chit[0].end = sc + 17;
     chit[0].anticodon = (sc[7] << 4) + (sc[8] << 2) + sc[9];	
     if (bp[sc[-1]][sc[17]])
      { chit[1].pos = sc-1;
        chit[1].stem = 6;
        chit[1].loop = 7;
        chit[1].looppos = sc + 5;
        chit[1].arm = 19;
        chit[1].end = sc + 18;
        chit[1].anticodon = chit[0].anticodon;
        nch = 2; }
     else nch = 1;

   /* 6 base cloop */
   /* exclude cstem < 4 Watson-Crick basepairs or equivalent */
   /* exclude cloop = RRnnNN */
   /* exclude cloop = NNnnYY */

     CLOOP6:
     if (cloop7) goto CLOOPE;
     if ((cm & 0xf00) >= 0x800) 
      { if (!YI[sc[6]])
	 if (!YI[sc[5]])
          goto CLOOP8;
	if (!RI[sc[9]])
         if (!RI[sc[10]])
	  goto CLOOP8;
        se = sc + 20;
	sg = sc;
	while (sg < se)
	 { sf = sg + 5;
	   while (sf < (sg + 11)) 
	    { if (*sf == *sg)
               if (sf[1] == sg[1])
                if (sf[2] == sg[2])
                 if (sf[3] == sg[3])
                  if (sf[4] == sg[4])
	           { sb = sg + 5;
		     s = sf + 5;
		     i = 0;
		     while (sb < sf) 
		      if (*sb++ != *s++) 
		       if (++i > 1) goto NXSEG6;
                     goto CLOOPE; }
	      NXSEG6:
	      sf++; }
	   sg++; }
	chit[nch].pos = sc;
        chit[nch].stem = 5;
        chit[nch].loop = 6;
        chit[nch].looppos = sc + 5;
        chit[nch].arm = 16;
        chit[nch].end = sc + 16;
        chit[nch++].anticodon = 0;
        if (athresh > 10) athresh = 10;
        if (bp[sc[-1]][sc[16]])
         { chit[nch].pos = sc-1;
           chit[nch].stem = 6;
           chit[nch].loop = 6;
           chit[nch].looppos = sc + 5;
           chit[nch].arm = 18;
           chit[nch].end = sc + 17;
           chit[nch++].anticodon = 0; }}

   /* 8 base cloop */
   /* exclude cstem < 4 Watson-Crick basepairs or equivalent */
   /* exclude cloop = RRnnnnNN */
   /* exclude cloop = NNnnnnYY */

     CLOOP8:
     if ((cm & 0xf) >= 0x8)
      { if (!YI[sc[5]])
	 if (!YI[sc[6]])
          goto CLOOPE;
	if (!RI[sc[12]])
         if (!RI[sc[11]])
	  goto CLOOPE;
        se = sc + 20;
	sg = sc;
	while (sg < se)
	 { sf = sg + 5;
	   while (sf < (sg + 11)) 
	    { if (*sf == *sg)
               if (sf[1] == sg[1])
                if (sf[2] == sg[2])
                 if (sf[3] == sg[3])
                  if (sf[4] == sg[4])
	           { sb = sg + 5;
		     s = sf + 5;
		     i = 0;
		     while (sb < sf) 
		      if (*sb++ != *s++) 
		       if (++i > 1) goto NXSEG8;
                     goto CLOOPE; }
	      NXSEG8:
	      sf++; }
	   sg++; }
        chit[nch].pos = sc;
        chit[nch].stem = 5;
        chit[nch].loop = 8;
        chit[nch].looppos = sc + 5;
        chit[nch].arm = 18;
        chit[nch].end = sc + 18;
        chit[nch++].anticodon = 0;
        if (athresh > 10) athresh = 10;
        if (bp[sc[-1]][sc[18]])
         { chit[nch].pos = sc-1;
           chit[nch].stem = 6;
           chit[nch].loop = 8;
           chit[nch].looppos = sc + 5;
           chit[nch].arm = 20;
           chit[nch].end = sc + 19;
           chit[nch++].anticodon = 0; }}

  /* calculate carm energy */

     CLOOPE:
     if (nch < 1) continue; 
     for (nc = 0; nc < nch; nc++)
           { s1 = chit[nc].pos;
             cstem = chit[nc].stem;
             cloop = chit[nc].loop;
             s4 = s1 + cstem;
             s2 = s4 + cloop;
	     energy = (cloop == 7)?0.0:-4.0;
             energy += cY[*s4] + cT[s4[1]] + cR[s2[-2]] + cA[s2[-1]];
             if (s4[1] == Cytosine)
              if (*s4 == Cytosine)
               energy += 2.6;
              else
               if (*s4 == Thymine)
                if (s2[-1] == Adenine)
                 energy += 2.6;
             s2 += cstem;
           stem_energy = bem[*s1][*--s2];
           k = neighbour_map[*s1][*s2];
           stem_energy += neighbour_em[k][s1[1]][s2[-1]];
           bondtype = btmap[*s1][*s2];
           if (bp[*s1][*s2])
            { if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
              else stem_energy += send_em[*s2][*s1]; }
           else
            { if (assymst[*s2][*s1]) stem_energy += mtTERMSTAB;
              else stem_energy += send_em[s2[-1]][s1[1]]; }
           while (++s1 < s4)
          { if (!wcbp[*s1][*--s2])
             { if (!wcbp[s1[-1]][s2[1]])
    { for (j = 0; j < mtNTM; j++)
       if (*s1 == tandemid[j][1])
        if (*s2 == tandemid[j][3])
                       if (s1[-1] == tandemid[j][0])
                        if (s2[1] == tandemid[j][2])
           { stem_energy += tandem_em[j];
             break; }
      if (s1 < (s4-1))
       if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB; }
               k = neighbour_map[*s1][*s2];
               stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                     neighbour_em[k][s1[1]][s2[-1]]); }
           bondtype += btmap[*s1][*s2];
             stem_energy += bem[*s1][*s2]; }
      if (!bp[*--s1][*s2])
       { s1--;
         s2++; }
      if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
             else stem_energy += send_em[*s1][*s2];
             cgcc = bondtype & 0xf;
             if (cgcc <= 0)
              { catc = (bondtype & 0xf0) >> 4;
                if (catc < cstem) energy -= mtGCPENALTY; }
             if (cstem == 6) energy += 1.0;
      chit[nc].bondtype = bondtype;
      chit[nc].stem_energy = stem_energy;
      chit[nc].energy = energy + stem_energy; }

  /* find tarms */

  nth = 0;
  slm = sc + 61; 
  sle = sc + 57;
  sb = sc + 21;
  sg = sc + 16;
  sge = sg + 30;
  slb = sg + 32;
  tem[0] = At[*slm];
  tem[1] = Ct[*slm];
  tem[2] = Gt[*slm];
  tem[3] = Tt[*slm];
  while (--slm > sle)
   { tem[0] = (tem[0] << 4) | At[*slm];
     tem[1] = (tem[1] << 4) | Ct[*slm];
     tem[2] = (tem[2] << 4) | Gt[*slm];
     tem[3] = (tem[3] << 4) | Tt[*slm]; }
  while (slm >= sb)
   { tem[0] = ((tem[0] << 4) | At[*slm]) & 0xfffff;
     tem[1] = ((tem[1] << 4) | Ct[*slm]) & 0xfffff;
     tem[2] = ((tem[2] << 4) | Gt[*slm]) & 0xfffff;
     tem[3] = ((tem[3] << 4) | Tt[*slm]) & 0xfffff;
     sf = slm + 3;
     if (sf > sge) sf = sge;
     apos2 = slm + 5;
     si = sg;
     s = si + 4;
     r = tem[*si];
     while (++si < s) r = (r >> 4) + tem[*si];
     while (si <= sf)
      { if (si < slm)
          r = (r >> 4) + tem[*si++];
        else
         { si++;
           r = r >> 4; }
     q = r & 0xf;
     if (slm > slb)
      { if (q < 5) continue;
        tloop = (int)(slm - si); }
     else
      { if (q < 2) continue;
        if (q < 3)
         { if (!wcbp[si[-5]][apos2[-1]]) continue;
           if (!wcbp[si[-4]][apos2[-2]]) continue;
           tloop = (int)(slm - si);
	   if (tloop > 5) continue; }
        else 
         { tloop = (int)(slm - si);
           if (q < 4)
            if (!bp[si[-4]][apos2[-2]])
             if (!bp[si[-2]][apos2[-4]])
              { if (tloop < 4) continue;
                if (si[-1] != Guanine) continue;
                if (*si != Thymine) continue;
                if (si[1] != Thymine) continue; }}}
     if (tloop < 7)
      { if (tloop < 2)
         if (tloop <= 0)
          { if (tloop <= -2)
             { if (!wcbp[si[-5]][apos2[-1]]) continue;
               if (!wcbp[si[-4]][apos2[-2]]) continue;
               tstem = 2;
               tloop += 6; }
            else
             if (bp[si[-3]][apos2[-3]])
              { tstem = 3;
                tloop += 4; }
             else
              { if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6; }}
         else
          { if (bp[si[-2]][apos2[-4]])
             { tstem = 4;
               tloop += 2; }
            else
             if (bp[si[-3]][apos2[-3]])
             { tstem = 3;
               tloop += 4; }
             else
              { if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6; }}
        else
         { if (bp[si[-1]][apos2[-5]])
            { if (q != 4) tstem = 5;
              else
               { if (bp[si[-2]][apos2[-4]]) tstem = 5;
                 else
                  { k = GI[si[-3]] + TI[si[-2]] + TI[si[-1]] + CI[*si];
                    if (k >= 2)
                     { tstem = 3;
                       tloop += 4; }
                    else tstem = 5; }}}
           else
            { if (bp[si[-2]][apos2[-4]])
               { tstem = 4;
                 tloop += 2; }
              else
               if (bp[si[-3]][apos2[-3]])
                { tstem = 3;
                  tloop += 4; }
               else
                { if (!wcbp[si[-5]][apos2[-1]]) continue;
                  if (!wcbp[si[-4]][apos2[-2]]) continue;
                  tstem = 2;
                  tloop += 6; }
	    }}
        if (tloop < 3)
         if (tstem > 3)
          { tstem--;
            tloop += 2; }}
     else
      { if (!bp[si[-1]][apos2[-5]])
         if (!bp[si[-2]][apos2[-4]])
          { tstem = 3;
            tloop += 4; }
         else
          { tstem = 4;
            tloop += 2; }
        else tstem = 5; }
     if (tloop > 17)
      if (tstem < 5)
       continue;

  /* calculate tarm energy */

     s1 = si - 5;
     tpos = s1;
     s4 = s1 + tstem;
     s2 = apos2;
     if (tt[*s1][*--s2])
      { energy = mtTSTTSTAB;
        if (tt[*++s1][*--s2])
	 { energy += mtTSTTSTAB;
	   bondtype = btmap[*s1++][*s2--]; }
	else bondtype = 0; }
     else
      { energy = 0.0;
	bondtype = 0; }

  /* calculate tstem energy */

     stem_energy = bem[*s1][*s2];
     k = neighbour_map[*s1][*s2];
     stem_energy += neighbour_em[k][s1[1]][s2[-1]];
     bondtype += btmap[*s1][*s2];
     while (++s1 < s4)
      { if (!wcbp[*s1][*--s2])
         { if (!wcbp[s1[-1]][s2[1]])
            { for (j = 0; j < mtNTM; j++)
              if (*s1 == tandemid[j][1])
               if (*s2 == tandemid[j][3])
                if (s1[-1] == tandemid[j][0])
                 if (s2[1] == tandemid[j][2])
                  { stem_energy += tandem_em[j];
                    break; }
              if (s1 < (s4-1))
               if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB; }
           k = neighbour_map[*s1][*s2];
           stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                     neighbour_em[k][s1[1]][s2[-1]]); }
        bondtype += btmap[*s1][*s2];
        stem_energy += bem[*s1][*s2]; }
      s1--;
      if (tloop < 4) stem_energy += ssend_em[*s1][*s2];
      else
       if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
       else stem_energy += send_em[*s1][*s2];

  /* compile possible tarms */

      energy += (stem_energy - mtBONDSTAB*(double)(5-tstem));
      if (energy >= tarmthresh)
       { thit[nth].pos = tpos;
         s1 = tpos + tstem;
	 s2 = apos2 - tstem;
         thit[nth].energy = energy - loop_stab[tloop] +
                  tG[s1[-1]] + tT[*s1] + tT[s1[1]] + tC[s1[2]];
         thit[nth].stem_energy = stem_energy;
         thit[nth].bondtype = bondtype;
         thit[nth].stem = tstem;
         thit[nth].loop = tloop;
	 thit[nth].end = tpos + 2*tstem + tloop;
         if (++nth >= mtNTH)
          { fprintf(stderr,"Too many mt-tstem hits\n");
            break; }
         if (tstem > 2)
	  if (tloop < 10)
           if (gt[s1[-1]][*s2])
           { thit[nth].pos = tpos;
             thit[nth].energy = energy - mtBONDSTAB - mtGTBOND -
		                loop_stab[tloop+2] +
                                tG[s1[-2]] + tT[s1[-1]] + tT[*s1] + tC[s1[1]];
             thit[nth].stem_energy = stem_energy - mtGTBOND;
             thit[nth].bondtype = bondtype - 0x100;
             thit[nth].stem = tstem - 1;
             thit[nth].loop = tloop + 2;
	     thit[nth].end = thit[nth-1].end;
             if (++nth >= mtNTH)
              { fprintf(stderr,"Too many mt-tstem hits\n");
                break; }
             if (tstem > 3)
	      if (tloop < 8)
	       if (gt[s1[-2]][s2[1]])
               { thit[nth].pos = tpos;
                 thit[nth].energy = energy - 2.0*mtBONDSTAB - 2.0*mtGTBOND -
		                loop_stab[tloop+4] + tG[s1[-3]] +
				tT[s1[-2]] + tT[s1[-1]] + tC[*s1];
                 thit[nth].stem_energy = stem_energy - 2.0*mtGTBOND;
                 thit[nth].bondtype = bondtype - 0x200;
                 thit[nth].stem = tstem - 2;
                 thit[nth].loop = tloop + 4;
	         thit[nth].end = thit[nth-1].end;
                 if (++nth >= mtNTH)
                  { fprintf(stderr,"Too many mt-tstem hits\n");
                    break; }}}
         if (tstem < 5)
          { if (tloop < 11) continue;
            if (tloop > 16) continue;
            if (!wcbp[s1[1]][s2[-2]]) continue;
            bondtype += btmap[*s1][s2[-1]] + btmap[s1[1]][s2[-2]];
            tstem += 2;
            tloop -= 4; }
         else
          { if (tloop < 9) continue;
            if (wcbp[*s1][s2[-1]])
             { if (tloop > 14) continue;
               tstem++;
               tloop -= 2;
               bondtype += btmap[*s1][s2[-1]]; }
            else
             { if (tloop < 11) continue;
               if (tloop > 16) continue;
               if (!wcbp[s1[1]][s2[-2]]) continue;
               bondtype += btmap[*s1][s2[-1]] + btmap[s1[1]][s2[-2]];
               tstem += 2;
               tloop -= 4; }}
         thit[nth].pos = tpos;
         s1 = tpos + tstem;
         thit[nth].energy = energy - loop_stab[tloop] +
                  tG[s1[-1]] + tT[*s1] + tT[s1[1]] + tC[s1[2]];
         thit[nth].stem_energy = stem_energy;
         thit[nth].bondtype = bondtype;
         thit[nth].stem = tstem;
         thit[nth].loop = tloop;
         thit[nth].end = thit[nth-1].end;
         if (++nth >= mtNTH)
          { fprintf(stderr,"Too many mt-tstem hits\n");
            break; }
         if (tloop < 9) continue;
         if (!wcbp[*s1][apos2[-tstem-1]]) continue;
         if (++tstem > 7) continue;
	 if (tloop > 14) continue;
         tloop -= 2;
         thit[nth].pos = tpos;
         s1 = tpos + tstem;
         thit[nth].energy = energy - loop_stab[tloop] +
                  tG[s1[-1]] + tT[*s1] + tT[s1[1]] + tC[s1[2]];
         thit[nth].stem_energy = stem_energy;
         thit[nth].bondtype = bondtype;
         thit[nth].stem = tstem;
         thit[nth].loop = tloop;
         thit[nth].end = thit[nth-1].end;
         if (++nth >= mtNTH)
          { fprintf(stderr,"Too many mt-tstem hits\n");
            break; }}}
     slm--; }

  /* find darms */

  ndh = 0;
  sle = sc - 4;
  slb = sc - 8;
  slm = sc - 1;
  tem[0] = dAt[*slm];
  tem[1] = dCt[*slm];
  tem[2] = dGt[*slm];
  tem[3] = dTt[*slm];
  while (--slm > sle)
   { tem[0] = (tem[0] << 4) | dAt[*slm];
     tem[1] = (tem[1] << 4) | dCt[*slm];
     tem[2] = (tem[2] << 4) | dGt[*slm];
     tem[3] = (tem[3] << 4) | dTt[*slm]; }
  slm1 = slm;
  while (slm > slb)
   { tem[0] = ((tem[0] << 4) | dAt[*slm]) & 0xffff;
     tem[1] = ((tem[1] << 4) | dCt[*slm]) & 0xffff;
     tem[2] = ((tem[2] << 4) | dGt[*slm]) & 0xffff;
     tem[3] = ((tem[3] << 4) | dTt[*slm]) & 0xffff;
     slm--;
     si = slm - 18;
     s = si + 3;
     r = tem[*si];
     while (++si < s) r = (r >> 4) + tem[*si];
     while (si <= slm1)
      { if (si < slm) r = (r >> 4) + tem[*si++];
        else
         { r = r >> 4;
	   si++; }
        if ((q = (r & 0xf)) < 6)
         { q += (unsigned int)(TI[si[-6]] + RI[si[-5]]);
           if (q < 6) continue; }

  /* calculate darm energy */

        s1 = si - 4;
        dhit[ndh].pos = s1;
        energy = dT[s1[-2]] + dA[s1[-1]];
        dloop = (int)(slm1 - si);
        if (dloop > 2)
         if (bp[si[-1]][*slm1])
          { dstem = 4;
	    goto EC; }
        if (dloop > 0)
         if ((ggstembp[si[-2]][slm[2]]) || (gabp[si[-1]][*slm1]))
            { dstem = 3;
              dloop += 2;
              energy += mtNOBOND;
	      goto EC; }
	if (!wcbp[si[-3]][slm[3]]) continue;
	if (!gc[si[-4]][slm[4]]) continue;
        dstem = 2;
        dloop += 4;
	if (dloop > 5) energy += mtNOBOND;
        energy += mtNOBOND;

        EC:
        s2 = slm + 4;
        s4 = s1 + dstem;
        if (!wcbp[s1[1]][s2[-1]])
         if (stemterm[s1[1]][s2[-1]]) energy -= 1.0;
         else 
          if (bp[s1[1]][s2[-1]]) energy -= 1.5;
          else energy -= 2.0;

  /* calculate dstem energy */

        stem_energy = bem[*s1][*s2];
        k = neighbour_map[*s1][*s2];
        stem_energy += neighbour_em[k][s1[1]][s2[-1]];
        bondtype = btmap[*s1][*s2];
        if (bp[*s1][*s2])
         { if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
           else stem_energy += send_em[*s2][*s1];
           s1++;
           s2--; }
        else
         { s1++;
           s2--;
           if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
           else stem_energy += send_em[*s2][*s1]; }
        stem_energy += bem[*s1][*s2];
        k = neighbour_map[*s1][*s2];
        stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                       neighbour_em[k][s1[1]][s2[-1]]);
        bondtype += btmap[*s1][*s2];
        while (++s1 < s4)
         { if (!wcbp[*s1][*--s2])
            { if (!wcbp[s1[-1]][s2[1]])
               { for (j = 0; j < mtNTM; j++)
                  if (*s1 == tandemid[j][1])
                   if (*s2 == tandemid[j][3])
                    if (s1[-1] == tandemid[j][0])
                     if (s2[1] == tandemid[j][2])
                      { stem_energy += tandem_em[j];
                        break; }
                 if (s1 < (s4-1))
                  if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB; }
              k = neighbour_map[*s1][*s2];
              stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                              neighbour_em[k][s1[1]][s2[-1]]); }
           bondtype += btmap[*s1][*s2];
           stem_energy += bem[*s1][*s2]; }
        if (!bp[*--s1][*s2])
         { s1--;
           s2++; }
        if (dloop < 4) stem_energy += ssend_em[*s1][*s2];
        else
         if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
         else stem_energy += send_em[*s1][*s2];

  /* compile possible darms */

              energy += stem_energy;
              dhit[ndh].energy = energy;
              dhit[ndh].stem_energy = stem_energy;
              dhit[ndh].bondtype = bondtype;
              dhit[ndh].stem = dstem;
              dhit[ndh].loop = dloop;
              if (++ndh >= mtND)
               { fprintf(stderr,"Too many mt-dstem hits\n");
                 break; }
	      if (dstem == 4)
	       { if (dloop >= 6)
	          if (bondtype < 0x1000)
               { s1 = si - 5;
	         s2 = slm + 5;
                 if (bp[*s1][*s2])
                  { dhit[ndh].pos = s1;
		    e = 0.5 + bem[*s1][*s2];
                    dhit[ndh].energy = energy + e;
		    if (wcbp[*s1][*s2])
		     dhit[ndh].energy += (dT[s1[-2]] + dA[s1[-1]] -
			                  dT[s1[-1]] - dA[*s1]);
                    dhit[ndh].stem_energy = stem_energy + e;
                    dhit[ndh].bondtype = bondtype + btmap[*s1][*s2];
                    dhit[ndh].stem = 5;
                    dhit[ndh].loop = dloop;
                    if (++ndh >= mtND)
                     { fprintf(stderr,"Too many mt-dstem hits\n");
                       break; }}}}
              else
               if (dloop >= 6)
               { s1 = si - 1;
	         s2 = slm1;
                 if (stemterm[*s1][*s2])
                  { dhit[ndh].pos = si - 4;
                    dhit[ndh].energy = energy;
                    dhit[ndh].stem_energy = stem_energy;
                    dhit[ndh].bondtype = bondtype;
                    dhit[ndh].stem = 4;
                    dhit[ndh].loop = dloop - 2;
                    if (++ndh >= mtND)
                     { fprintf(stderr,"Too many mt-dstem hits\n");
                       break; }}}
              if (dloop >= 4) continue;
              s1 = si - 4 + dstem - 1;
              s2 = s1 + dloop + 1;
              if (bp[*s1][*s2]) continue;
              dhit[ndh].pos = si - 4;
              dhit[ndh].energy = energy + 0.001;
              dhit[ndh].stem_energy = stem_energy;
              dhit[ndh].bondtype = bondtype;
              dhit[ndh].stem = dstem - 1;
              dhit[ndh].loop = dloop + 2;
              if (++ndh >= mtND)
               { fprintf(stderr,"Too many mt-dstem hits\n");
                 break; } }
     slm1--; }

  /* build darm exclusion map */
  /* 5' astems further from carm than */
  /* mt_DRLmaxlength must match a darm */

  for (i = 3; i <= 30; i++) dposmap[i] = 0;
  sf = sc - mt_DRLmaxlength - 1;
  sld = sf;
  if (ndh > 0)
   { s = dhit[0].pos;
     for (nd = 0; nd < ndh; nd++)
      { se = dhit[nd].pos;
        if (se < s) s = se;
	i = (int)(sc - se);
	if (dposmap[++i] < 1) dposmap[i] = 1;
	dposmap[++i] = 2;
	if (dposmap[++i] < 1) dposmap[i] = 1; }
     s -= 4;
     if (s < sf) sf = s; }
	
  /* build tarm exclusion map */
  /* 3' astems further from carm than */
  /* mt_TVRLmaxlength must match a tarm */

  for (i = 17; i <= 62; i++) tendmap[i] = 0;
  s2 = sc + mt_TVRLmaxlength + 17;
  sle = s2;
  if (nth > 0)
   { s = thit[0].end;
     for (nt = 0; nt < nth; nt++)
      { se = thit[nt].end;
	    if (se > s) s = se;
	    i = (int)(se - sc);
	    bondtype = thit[nt].bondtype;
	    if (tendmap[i])
         { if (bondtype < tendmap[i]) tendmap[i] = bondtype; }
	    else tendmap[i] = bondtype; }
     if (s > s2) s2 = s; }

  /* find astems in 3 categories: */
  /* high energy astems close to carm */
  /* high energy astems matching a high energy tarm far from carm */
  /* low energy astem matching a darm and tarm */

     nah = 0;
     sa = sc - 3;
     sg = sf - 6;
     sb = sc + 17;
     se = s2 + 6;
     tem[0] = aAt[*se];
     tem[1] = aCt[*se];
     tem[2] = aGt[*se];
     tem[3] = aTt[*se];
     while (--se > s2)
      { tem[0] = (tem[0] << 4) | aAt[*se];
        tem[1] = (tem[1] << 4) | aCt[*se];
        tem[2] = (tem[2] << 4) | aGt[*se];
        tem[3] = (tem[3] << 4) | aTt[*se]; }
     ti = (int)(se - sc);
     while (se >= sb)
      { tem[0] = ((tem[0] << 4) | aAt[*se]) & 0xfffffff;
        tem[1] = ((tem[1] << 4) | aCt[*se]) & 0xfffffff;
        tem[2] = ((tem[2] << 4) | aGt[*se]) & 0xfffffff;
        tem[3] = ((tem[3] << 4) | aTt[*se]) & 0xfffffff;
	if (tendmap[ti])
	 { nti = (tendmap[ti] < 0x2000)?1:0; }
	else
	 { if (se > sle) goto ANX;
	   nti = -1; }
        si = sg;
        r = tem[*si];
        while (++si < sf) r = (r >> 4) + tem[*si];
		di = (int)(sc - si);
        while (si < sa)
         { r = (r >> 4) + tem[*si++];
	   if (dposmap[--di])
            { if (nti <= 0)
               { if (nti < 0)
                  if (dposmap[di] < 2) continue;
	         if ((av = (r & 0xf)) < athresh) continue; }}
	   else
	    { if (si < sld) continue;
	      if (nti < 0) continue;
              if ((av = (r & 0xf)) < athresh) continue; }
           if (nah >= mtNA)
            { fprintf(stderr,"Too many mt-astem hits\n");
              break; }

  /* predict astem length and calculate astem energy */

       s1 = si - 7;
       s2 = se + 6;
       if (bp[*s1][*s2])
        { astem = 7;
          energy = 0.0;
          ahit[nah].pos1 = s1;
          ahit[nah].pos2 = se; }
       else
        if (ggstemterm[*s1][*s2])
         { astem = 7;
           ahit[nah].pos1 = s1;
           ahit[nah].pos2 = se;
           energy = bem[*s1++][*s2--]; }
	else
         { energy = bem[*s1++][*s2--];
           if (bp[*s1][*s2])
            { astem = 6;
              ahit[nah].pos1 = s1;
              ahit[nah].pos2 = se; }
           else
            if (ggstemterm[*s1][*s2])
             { astem = 6;
               ahit[nah].pos1 = s1;
               ahit[nah].pos2 = se;
               energy += bem[*s1++][*s2--]; }
	    else
             { astem = 5;
               energy += bem[*s1++][*s2--];
               ahit[nah].pos1 = s1;
               ahit[nah].pos2 = se; }}
       ahit[nah].stem = astem;
       bondtype = btmap[*s1][*s2];
       energy += bem[*s1][*s2];
       k = neighbour_map[*s1][*s2];
       energy += neighbour_em[k][s1[1]][s2[-1]];
       energy += bem[*++s1][*--s2];
       k = neighbour_map[*s1][*s2];
       energy += (neighbour_em[k][s1[-1]][s2[1]] +
                 neighbour_em[k][s1[1]][s2[-1]]);
       bondtype += btmap[*s1][*s2];
       while (++s1 < si)
        { if (!wcbp[*s1][*--s2])
           { if (!wcbp[s1[-1]][s2[1]])
              { for (j = 0; j < mtNTM; j++)
                 if (*s1 == tandemid[j][1])
                  if (*s2 == tandemid[j][3])
                   if (s1[-1] == tandemid[j][0])
                    if (s2[1] == tandemid[j][2])
                     { energy += tandem_em[j];
                       break; }
                if (s1 < (si-1))
                 if (!bp[s1[1]][s2[-1]]) energy -= mt3MMSTAB; }
             k = neighbour_map[*s1][*s2];
             energy += (neighbour_em[k][s1[-1]][s2[1]] +
             neighbour_em[k][s1[1]][s2[-1]]); }
          bondtype += btmap[*s1][*s2];
          energy += bem[*s1][*s2]; }
       if (!bp[*--s1][*s2])
        if (!bp[*--s1][*++s2])
         if (!bp[*--s1][*++s2])
          if (!bp[*--s1][*++s2])
           goto NOST;
       if (assymst[s1[1]][s2[-1]]) energy += mtTERMSTAB;
       NOST:
       ahit[nah].energy = energy;
       ahit[nah].bondtype = bondtype;
       nah++; }
      ANX:
      se--;
      ti--; }
    if (nah <= 0) continue;

	
  /* build mttrna genes */
  /* cycle through astems first so that */
  /* GC content is only calculated once per astem */

  thresh = -INACTIVE;
  te.ps = NULL;
  for (na = 0; na < nah; na++)
      { apos2 = ahit[na].pos2;
        apos1 = ahit[na].pos1;
        astem = ahit[na].stem;
        aend1 = apos1 + astem;
        astem8 = (astem == 7)?(wcbp[apos1[-1]][apos2[7]]):0;
        asteme = 0;
        ea = ahit[na].energy;
        abondtype = ahit[na].bondtype;
	agcat = ((abondtype >> 4) + abondtype) & 0xf;

  /* GC content */

    s = apos1;
	aend2 = apos2 + astem;
	nbase = (int)(aend2 - apos1) + 1;
	igc = 0;
	while (s <= aend2)
         { k = *s++;
           if (k >= Cytosine) if (k <= Guanine) igc++; }
        gcv = 10.0*(double)igc/(double)nbase;
	if (gcv < 1.0)
         { if (gcv < 0.55) continue;
	   ea -= 0.5; }
	if (nbase > 60)
	 { if (gcv > 6.0) ea -= 2.0*(gcv - 6.0); }
	else
	 { if (gcv > 5.0) ea -= 2.0*(gcv - 5.0); }
	if (gcv > 6.6)
         { ea -= 6.0;
           if (gcv > 7.0) ea -= 6.0; }

  /* findout if inside a coding sequence */



        incds = 0;
        i = -1;
        while (++i < ncdsh)
         if (apos1 > cdshit[i].pos1)
          if (aend2 <= cdshit[i].pos2)
           { incds = 1;
             ea -= 2.0;  
             break; }


  /* cycle through carms that fall between astem */

        nc = -1;
        while (++nc < nch)
         { cpos = chit[nc].pos;
           dloop = (int)(cpos - aend1);
           if (dloop < 3) continue;
           if (dloop > 26) continue;
           cend = chit[nc].end;
           tloop = (int)(apos2 - cend);
           if (tloop < 5) continue;
	   cloop = chit[nc].loop;
           cstem = chit[nc].stem;
           clooppos = chit[nc].looppos;
	   cloopend = clooppos + cloop;
           carm = chit[nc].arm;
           anticodon = chit[nc].anticodon;
	   cbondtype = chit[nc].bondtype;
	   acbondtype = abondtype + cbondtype;
	   cgcat = ((cbondtype >> 4) + cbondtype) & 0xf;
           ec = ea + chit[nc].energy;

  /* astem,cstem stability (GC bond count) */

         if ((abondtype & 0xf) <= 0)
          if ((cbondtype & 0xf) <= 0)
           { ec -= mtGCPENALTYD;
             if (((cbondtype & 0xf0) >> 4) >= 5) ec += 0.5; }
		

  /* anticodon to astem discriminator base match */

         astem8d = 0;
         if (cloop == 7)
          { if (!mt_discrim[ds][anticodon][apos2[astem]])
             if (astem8)
              if (mt_discrim[ds][anticodon][apos2[8]]) astem8d = 1;
              else ec -= 3.0;
             else
              if (astem <= 6)
               { if (!mt_discrim[ds][anticodon][apos2[7]])
                  if (astem == 5)
                   { if (!mt_discrim[ds][anticodon][apos2[6]]) ec -= 3.0; }
                  else
                   ec -= 3.0; }
              else
               ec -= 3.0; }
		
		
  /* build TV-replacement loop mttrna genes */

        if (tloop <= mt_TVRLmaxlength)
         { if (!sw->tvloop) goto TVN;

  /* astem termination */
  /* (only need to calculate once per astem) */

    if (!asteme)
            { asteme = 1;
              s = aend1 - 1;
              se = apos2;
              while (!bp[*s][*se])
               { if (--s <= apos1)
		  { eas = 0.0;
		    goto NOST2; }
                 se++; }
              if (!aastemterm[s[1]][se[-1]]) eas = -0.5;
      else
       { eas = 0.0;
         while (se >= apos2)
                 { s++;
                   se--;
         if (aastemterm[*s][*se]) eas += 1.0;  }}}

  /* choose darm */

        NOST2:
           energy = 94.0 + ec + eas;
           nd = -1;
           ndi = -1;
           ed = -INACTIVE;
           while (++nd < ndh)
            { dpos = dhit[nd].pos;
              spacer1 = (int)(dpos - aend1);
              if (spacer1 != 2) continue;
              dl = dhit[nd].loop;
              dstem = dhit[nd].stem;
	      if (dstem > 4) continue;
              darm = 2*dstem + dl;
              spacer2 = (int)(cpos - dpos) - darm;

  /* astem,darm,cstem interspacing */

              if (spacer2 < 1) continue;
              e = dhit[nd].energy;
              if (spacer2 > 1)
               { if (spacer2 > 2) continue;
                 if (!stembp[*cpos][cend[-1]]) continue; 
		         if (tloop > 12) e -= 2.0;
	             if ((dhit[nd].bondtype & 0xf) < 1) 
		         if ((agcat + cgcat + 1) < (cstem + astem)) e -= 3.0; }
              else
               if (dl > 11)
	            { if (!RI[cpos[-1]]) e -= 2.0; }
               else
                { if (cpos[-1] == Cytosine) e -= 2.0; }

  /* small,large dloop, dstem R motif */

              if (dl < 3) e -= 2.0;
	          if (dl > 12) e -= 2.0;
              if (!RI[*dpos]) e -= 1.0;

  /* darm,tloop tertiary interaction */

              k = 0;
	          di = ((dl >= 12)?3:((dl >= 9)?2:1));
              tl = (tloop >= 14)?5:((dl >= 9)?((tloop >= 10)?4:3):3);
              if (!ggstackbp[dpos[dstem+di]][cend[tl]]) 
	           { if (tl > 3)
	              { if (!ggstackbp[dpos[dstem+di]][cend[tl-1]]) e -= 1.5;
	                else k++; }
                 else 
	              if (di > 1)
	               { if (!ggstackbp[dpos[dstem+di-1]][cend[tl]]) e -= 1.5;
	                 else k++; }
		          else
                    e -= 1.5; }   
              else k++;
              if (stemterm[dpos[dstem-1]][dpos[darm-dstem]])
               { e -= 0.5;
                 if (cend[2] == dpos[dstem-2])
                  { if (bp[cend[2]][dpos[darm-dstem+1]]) k++; }
                 else
                  { if (cend[2] == dpos[darm-dstem+1]) 
                     if (bp[cend[2]][dpos[dstem-2]]) k++; }}
              else
               { if (cend[2] == dpos[dstem-1])
                { if (!bp[cend[2]][dpos[darm-dstem]]) e -= 0.5;
                  else k++; }
               else
                { if (cend[2] != dpos[darm-dstem]) e -= 0.5;
                  else
                   if (!bp[cend[2]][dpos[dstem-1]]) e -= 0.5;
                   else k++; }}
              if (cend[1] == *dpos)
               { if (!stackbp[cend[1]][dpos[darm-1]]) e -= 0.5;
                 else k++; }
              else
               { if (cend[1] != dpos[darm-1]) e -= 0.5;
                 else
                  if (!bp[cend[1]][*dpos]) e -= 0.5;
                  else k++; }

  /* darm stability */

              dstemmotif = wcbp[dpos[1]][dpos[darm-2]];
              if (spacer2 == 2)
               if ((k < 3) || (dhit[nd].bondtype > 0x200) || (!dstemmotif))
                { if (abondtype >= 0x10000) e -= 2.0;
                  if (dstem > 3) e -= 1.0; 
		          e -= 0.5; }

  /* darm tertiary interactions */

              j = 0;
              b8 = dpos[-2];
              b9 = dpos[-1];
              if (!bp[b8][dpos[dstem]]) e -= 1.0;
              else if (wcbp[b8][dpos[dstem]]) j++;
              if (!bp[b8][dpos[darm-dstem-1]]) e-= 1.0;
              else if (wcbp[b8][dpos[darm-dstem-1]]) j++;
              if (!wcbp[dpos[2]][dpos[darm-3]])
               { if (!gastembp[b8][dpos[dstem]]) e -= 2.0;
                 else if (!gastembp[b8][dpos[darm-dstem-1]]) e -= 2.0;
	         if (!ggstembp[dpos[2]][dpos[darm-3]]) e -= 1.0; }
              else j++;
              if (!bp[b9][dpos[2]])
               { if (!bp[b9][dpos[darm-3]]) e -= 1.0;
                 else j++; }
              else j++;

/* more extensive tertiary interaction between darm,tloop */

	      if (dstemmotif)
               { if (k >= 3)
                  if (bp[dpos[2]][dpos[darm-3]])
                   { if (b8 != Thymine) e += 0.5;
                     if (dl > 3)
                      if (bp[dpos[dstem+2]][cend[tl+1]]) e += 0.7;
                      else
                       if (gabp[dpos[dstem+2]][cend[tl+1]]) e += 0.5;
		     if (tloop >= 6)
		      if (spacer2 < 2)
                       if (dl >= 3)
		        { di = (dl > 11)?2:1;
                          if (bp[dpos[dstem+di]][cend[tl]])
		           { if (chit[nc].stem_energy > -4.8) e += 0.5;
                             if (wcbp[dpos[dstem+di]][cend[tl]])
                              if (gcv > 1.2)
                               if (clooppos[1] == Thymine)
			        if (cbondtype < 0x200)
			         if ((cbondtype & 0xf) > 0)
                                  if (abondtype < 0x2000)
				   { e += 1.5;
				     if (dl > 3)
                                      if (wcbp[dpos[dstem+di+1]][cend[tl+1]])
                                       e += 1.0; }}}}
                 if (j >= 4) e += 0.25; }
              if (e > ed)
               { ed = e;
                 ndi = nd;
	         ti = k; }}
      if (ndi < 0) goto TVN;
      energy += ed;
      dpos = dhit[ndi].pos;
      dstem = dhit[ndi].stem;
      dl = dhit[ndi].loop;
      darm = 2*dstem + dl;
      dbondtype = dhit[ndi].bondtype;
      spacer2 = (int)(cpos - dpos) - darm;
      spacer1 = (int)(dpos - aend1);
      b8 = *aend1;
      b9 = aend1[1];

  /* false positive suppression */

      if (dloop < 15) energy -= 2.0;
      if (cstem > 5) energy -= 1.0;
      if (tloop < 6) energy -= 1.0;
      if (tloop > 12)
       { energy -= 1.0;
         if (agcat < 6) energy -= 2.0;
         if (tloop > 15) energy -= 2.5; }
      if (!stackbp[*dpos][dpos[darm-1]]) energy -= 1.0;
      if (dstem < 4)
       if (gcv > 1.2)
        if ((dbondtype & 0xf0f) == 0) energy -= 1.5;
     if (b8 != Thymine)
      { if (dl < 4)
	     if (abondtype > 0x10000)
           energy -= 1.5;
        if (b8 == Adenine) 
         if (YI[cloopend[-2]])
          energy -= 1.0; }  
      if (dl > 10)
       { if (tloop < 7) energy -= 2.0;
         if (spacer2 > 1) energy -= 2.0;
	 if (dhit[ndi].stem_energy < -3.4) energy -= 2.0; }
      if (gcv < 2.0)
       if (dbondtype > 0x10000) energy -= 2.0;
      if ((cbondtype & 0xf) < 1)
       if (abondtype > 0x100)
        { if (cgcat < 4)
           energy -= 1.5;
	  if (!wcbp[dpos[2]][dpos[darm-3]]) energy -= 1.0; }
      if (b8 != Thymine)
       if ((clooppos[1] != Thymine) ||
           (*clooppos != Cytosine))
        if (dl > 3)
         if (dbondtype > 0x10000)
           energy -= 1.0;
      if (!RI[cend[1]])
       if (b9 != Guanine) energy -= 1.0;
       else energy -= 0.5;
      if (b9 == Guanine)
       { if (!RI[*cend]) energy -= 1.0;
         if (spacer2 != 1) energy -= 3.0;
         else
          { tl = (tloop >= 14)?5:((dl >= 9)?((tloop >= 7)?4:3):3);
            s = dpos + dstem;
            if (!wcbp[s[1]][cend[tl]])
	     { energy -= 2.5;
               if (dl >= 5)
                if (chit[nc].energy > 2.0)
                 if (wcbp[s[2]][cend[tl]])
                  if (wcbp[s[3]][cend[tl+1]])
                   energy += 6.0; }
	    else
	     if (b8 == Thymine)
             if (dl >= 5)
              if (chit[nc].energy > 2.0)
               if (wcbp[s[2]][cend[tl+1]])
                   energy += 3.5; }}
      else if (b9 != Adenine) energy -= 3.0;
      if (b8 != Thymine)
       if (b8 == Guanine)
        { if (!RI[dpos[dstem]]) energy -= 1.0;
          else
           if (RI[dpos[darm-dstem-1]]) energy += 2.0; }
       else energy -= 1.0;

   /* carm termination */

       if (assymst[cend[-1]][*cpos]) energy += 1.0;

   /* CTnnnAA cloop motif */

       energy += CX7[*clooppos] + AX7[cloopend[-2]];
       if (clooppos[1] == Cytosine) energy -= 2.0;

   /* NNnnnAA cloop motif */

       if (cloopend[-2] == Adenine)
        if (cloopend[-1] == Adenine)
             if (spacer1 == 2)
              if (dbondtype < 0x1000)
               { if (abondtype < 0x100) energy += 1.0;
                 else
		  if (cbondtype < 0x100) energy += 1.0; }

  /* global stem damage level */

           bondtype = acbondtype + dbondtype;
           i = (int)((bondtype >> 16) & 0xf);
           j = (int)((bondtype >> 12) & 0xf);
           k = (int)((bondtype >> 8) & 0xf);
	   if (k > 0)
	    if (i > 0)
             { k += (i + j);
               if (k > 5) energy -= 1.0*(double)(k - 5); }

  /* global stem stability (GC bond count) */

      gcc = bondtype & 0xf;
      if (gcc < 2) 
       { if (ti >= 2)
          { if (cbondtype < 0x100)
             if ((cbondtype & 0xf) > 0) goto NGCC1;
            if (ti >= 3)
             if (cgcat >= 4)
              { if ((cbondtype & 0xf) > 0) goto NGCC1;
                if (cbondtype < 0x100) goto NGCC2; }}
	  energy -= (double)(3 - gcc);
	 NGCC2:
         if (gcc < 1)
          { if (agcat < 5) energy -= 2.0;
	    if (bondtype > 0x10000) energy -= 1.5; }}
	 NGCC1:


  /* global stability */
  /* (stem stability,dloop-tloop tertiary interaction,dloop size) */

	   if (abondtype > 0x1000)
            if (ti < 3)
	     { if (chit[nc].stem_energy < -6.0)
                energy -= 1.5;
               if (dl > 9)
		if (((dbondtype + cbondtype) & 0xf) < 1)
                 energy -= 1.0; }

  /* tloop,dloop tertiary interaction */
  /* (alternative dloop position) */

	   if (bondtype < 0x1000)
            if (b8 == Thymine)
             if (RI[b9])
              if (dl > 4)
               if (!bp[cend[3]][dpos[dstem+1]])
                if (bp[cend[3]][dpos[dstem+2]])
                 energy += 0.5;


  /* "near perfect" TV-loop mttRNA: */
  /* darm-tloop tertiary interaction,low global stem damage, */
  /* TR motif at b8-9, good astem,darm,carm interspacing */

           if (ti >= 2)
            if (agcat >= 6)
             if (cbondtype < 0x100)
              if (dbondtype < 0x100)
               if (RI[b9])
               if (b8 == Thymine)
                if ((abondtype & 0xf) > 0)
                 if ((dbondtype & 0xf) > 0)
                  if (spacer1 == 2)
                   if (spacer2 == 1)
                    energy += 1.5;

  /* find exceptions */

       if (energy < dthresh)
        { if (!mtxdetect) goto TVN;
          if (incds) goto TVN;
	      if (energy < (thresh - 7.0)) goto TVN;
	      if (energy < (dthresh - 7.0)) goto TVN; 
	      if (nbase > 68) goto TVN;
          if (abondtype > 0x20100) goto TVN;
	      if (dl > 9) 
               { if (dl > 10) goto TVN;
                 if (dstem < 4) goto TVN;
                 if (dbondtype > 0x100) goto TVN; }
	      if (dstem > 4) goto TVN;
          if (b9 != Adenine) 
               { if (b9 != Guanine) goto TVN;
		  if (cbondtype > 0x100) goto TVN;
		  if (dbondtype > 0x200) goto TVN; }
	      if (cloop != 7) goto TVN;
          if (YI[cloopend[-2]]) goto TVN;
          if (b8 == Thymine)
                { if (apos2[-1] == Thymine)
                   if (apos2[-2] == Thymine)
                    if (tloop < 8)
                     if (tt[aend1[-1]][*apos2])
                      if (wcbp[dpos[2]][dpos[darm-3]])
                       if (((dbondtype + cbondtype) & 0xf) > 0)
			energy += 3.0; }
               else
                if (b8 == Adenine)
                 { if (apos2[-1] == Adenine)
                    if (apos2[-2] == Adenine)
		     { if (assymat[aend1[-1]][*apos2])
                        if (assymat[apos2[1]][aend1[-2]]) energy += 2.0;
		       if (agcat >= 5)
		        if (cgcat >= 4)
                         if (dbondtype < 0x100)
		          if (at[aend1[-1]][*apos2])
                           if (at[apos2[1]][aend1[-2]])
                            energy += 1.0; }
                   if (ti >= 3)
                    if (cgcat >= 4)
                     if (agcat >= 4)
                      if ((cbondtype & 0xf) > 0)
                       if ((abondtype & 0xf) > 1)
                      if (dbondtype < 0x200)
                       if (wcbp[dpos[1]][dpos[darm-2]])
                        if (clooppos[1] == Thymine)
                         if (YI[*clooppos])
                          if (RI[cloopend[-2]])
                           if (RI[cloopend[-1]])
                            energy += 5.0; }
	      if (bondtype < 0x100)
               { if (spacer2 == 1)
	          if (*clooppos == Cytosine)
		    if (clooppos[1] == Thymine)
                     if (cloopend[-2] == Adenine)
                      if (cloopend[-1] == Adenine)
                       energy += 2.0; }
              else
               { if (spacer2 == 1)
                     { if (b8 == Thymine)
		        if (dl > 3)
                         if (dbondtype < 0x200)
                         { if (cbondtype < 0x100)
                            { if (!bp[dpos[dstem+1]][cend[3]])
                               if (bp[dpos[dstem+1]][cend[4]])
                                energy += 2.0;
                              if (dbondtype < 0x100)
                               if (abondtype < 0x20000)
                                if (ti >= 2)
                                 if (dstem >= 3)
                                  if (tloop < 13)
                                   if ((cbondtype & 0xf) > 0)
                                    energy += 4.0; }}
                         else
                          if (dstem > 3)
                          if (dbondtype < 0x300)
			   { if (bondtype < 0x10000)
                              if (ti >= 3)
                               if ((acbondtype & 0xf) > 0)
                                if (wcbp[dpos[2]][dpos[darm-3]])
                                 energy += 4.0; }
                       if (tloop < 8)
		        { if (dbondtype < 0x200)
			   { if (cbondtype < 0x100)
                              if (ti >= 2)
                               {
                               if (wcbp[dpos[dstem+1]][cend[3]])
			        {
                                if (b8 == Thymine)
			         if (abondtype < 0x3000)
                                  energy += 5.0;
                                if (agcat >= 5)
			         if (gcv > 1.2)
				  if (RI[cloopend[-1]])
                                   energy += 7.0;
			        }
		               if (dbondtype < 0x100)
                                if (agcat >= 6)
                                 if (YI[*clooppos])
                                  if (clooppos[1] == Thymine)
                                   if (RI[cloopend[-2]])
                                    if (RI[cloopend[-1]])
                                     energy += 2.0; 
                                }
			     if (cbondtype < 0x300)
			      if (ti >= 3)
                               if (abondtype < 0x2000)
				if ((dbondtype & 0xf) > 0)
				 if ((acbondtype & 0xf) > 0)
                                  if (ahit[na].energy >= -7.0)
                                   if (dstem >= 4)
                                     energy += 3.0; }
			  if (dbondtype < 0x300)
                           if (cgcat >= 4)
                            if (abondtype < 0x2000)
                             if (ahit[na].energy >= -7.0)
                              if (cbondtype < 0x10000)
                               if ((cbondtype & 0xf) > 0)
                                if (cstem < 6)
			         if (ti >= 3)
			          energy += 4.0; }}
                  if (tloop > 8)
                   if (agcat >= 6)
                    if (cbondtype < 0x100)
                     if ((cbondtype & 0xf) > 0)
                      if (b8 == Thymine)
                       if (wcbp[dpos[dstem+1]][cend[3]])
                        if (wcbp[dpos[1]][dpos[darm-2]])
                         energy += 7.0; }
              if (dbondtype < 0x100)
               if (cgcat >= 4)
                if (agcat >= 5)
                        if (wcbp[dpos[1]][dpos[darm-2]])
                               if ((cbondtype & 0xf) > 0)
                               if ((abondtype & 0xf) > 0)
                               if ((dbondtype & 0xf) > 0)
				       energy += 0.5;
	      if (cbondtype < 0x100)
	       if (dbondtype < 0x200)
                if (agcat >= 5)
                 if (b8 == Thymine)
                  if (tloop < 8)
                   if (wcbp[dpos[1]][dpos[darm-2]])
                    if (wcbp[dpos[2]][dpos[darm-3]])
                     if ((cbondtype & 0xf) > 0)
                      if ((abondtype & 0xf) > 0)
                       if ((dbondtype & 0xf) > 0)
                        if (clooppos[1] == Thymine)
                         if (YI[*clooppos])
                          if (RI[cloopend[-2]])
                           energy += 3.0;
	      if (energy < dthresh) goto TVN;
              energy -= (0.9*(energy - dthresh) + 5.0); }

  /* remember fully formed TV-loop replacement mttRNA gene */
  /* if threshold reached */

         if (energy < thresh) goto TVN;
         te.energy = energy;
         thresh = energy;
         te.ps = apos1;
         te.dstem = dstem;
         te.dloop = dl;
         te.spacer1 = spacer1;
         te.spacer2 = spacer2;
         te.cstem = cstem;
         te.cloop = cloop;
         k = astem + spacer1 + darm + spacer2;
         te.anticodon = k + cstem + 2;
         te.nintron = 0;
         te.intron = 0;
         te.var = 0;
         te.varbp = 0;
         te.tstem = 0;
         te.tloop = tloop;
         te.nbase =  k + carm + tloop;
	 tastem = astem;
	 tastem8 = astem8;
	 tastem8d = astem8d;
	
  /* build D-replacement loop mttrna genes */

          TVN:
          if (tloop < 10) continue; }
          if (dloop > mt_DRLmaxlength) goto DN;
	      if (gcv < 1.2) goto DN;
          energy = 91.0 + ec;

  /* CCnnnAA cloop */

          if (clooppos[1] == Cytosine)
            { if (*clooppos != Cytosine) goto DN;
              if (cloopend[-2] != Adenine) goto DN;
              if (cloopend[-1] != Adenine) goto DN;
	      energy -= 1.0; }

  /* choose tarm */

           nt = -1;
           nti = -1;
           et = -INACTIVE;
           while (++nt < nth)
            { tl = thit[nt].loop;
              if (tl > 11) continue;
	      if (thit[nt].end != apos2) continue;
              tpos = thit[nt].pos;
              tstem = thit[nt].stem;

  /* var loop (3-7 bases long) */

              var = (int)(tpos - cend);
              if (var < 3) continue;
              e = thit[nt].energy;
              if (var > 5)
               { if (var > 7) continue;
		 if (tl < 7) continue;
		 e -= 1.0;
                 if ((dloop < 10) || (tstem < 4))
                  e -= 2.0*(double)(var - 5); }

  /* tloop RA or RG motif */

              s = tpos + tstem;
              k = 0;
	      n = 0;
              i = 0;
              while ((j = tloopa[tl][i++]) >= 0)
               if (s[j] == Adenine)
                { k = 1;
                  if (dloop >= 3)
                   if (tl > 3)
                    { b57 = s[j-1];
                      if (RI[b57] || (tl < 5))
                       { if (bp[b57][aend1[0]])
                          { e += 1.5;
                            n = 1;
                            break; }
                         if (bp[b57][aend1[1]])
                          { e += 1.5;
                            n = 1;
                            break; }
                         if (dloop > 10)
                          if (bp[b57][aend1[2]])
                           { e += 1.5;
                             n = 1;
                             break; }}}}
              if (!k)
               { i = 0;
                 while ((j = tloopa[tl][i++]) >= 0)
                  if (s[j] == Guanine)
                   if (RI[s[j-1]])
                    { k = 1;
                      break; }
                 if ( j < 0) e -= ((tl > 5)?2.0:1.0); }

  /* tertiary interaction between tloop and start of dloop */

              ti = (tl > 5)?1:((dloop > 5)?1:0);
              di = (dloop > 5)?2:1;
              if (stackbp[aend1[di]][s[ti]]) e += 1.0;

  /* tloop GTTC motif */

              i = (s[-1] == Guanine)?1:0;
              if (tl >= 5)
               { ti = i + TI[*s] + TI[s[1]] + CI[s[2]];
                 if (n)
                  if (!i)
                   if (TI[*s])
                    if (TI[s[1]])
                     if (AI[s[2]])
                      if (tl >= 7)
                       ti++;
                 if ((i > 0) || (ti >= 3)) e += (double)ti; }
              else
               { ti = i + TI[*s] + TI[s[1]];
                 if ((i > 0) || (ti >= 2)) e += (double)ti; }
              if (e > et)
               { et = e;
                 nti = nt;
                 tc = k; }}

        if (nti < 0) goto DN;
           energy += et;
           tpos = thit[nti].pos;
           tstem = thit[nti].stem;
           tl = thit[nti].loop;
           tbondtype = thit[nti].bondtype;
           var = (int)(tpos - cend);

  /* tertiary interaction between b48(=tpos[-1]) and dloop */

           b48 = tpos[-1];
           if (dloop <= 7)
            { if (YI[b48]) tc++;
              else energy -= 1.0; }
           else
            { i = 0;
              while ((j = dloopi[dloop][i++]) >= 0)
              if (assymagbp[b48][aend1[j]])
               { tc++;
                 break; }
              if (j < 0) energy -= 1.0; }

   /* large dloop, large tloop */

           if (dloop > 7)
	    { if (tl >= 6)
               if (tc < 2)
                energy -= 2.0;
	      if (tstem < 3)
                energy -= 1.0; }

   /* carm termination */

           s = cpos - 1;
           se = cend;
           if (cstem > 5)
            { s++;
              se--; }
           if (!stackbp[*s][*se]) energy -= 1.0;
           se = cpos - 3;
           if (!bp[cend[-1]][*cpos])
            { if (assymst[cend[-1]][*cpos])
               { if (dloop < 5) se++;
                 energy += 1.5; }
              else if (dloop < 13) se++; }
           else
            { if (cstem > 5) { if (dloop < 13) se++; }
              else if (dloop < 5) se++; }

   /* tertiary interaction between tloop and dloop near carm */

    s = tpos + tstem;
    if (tl >= 5)
     { ti = (tl >= 10)?4:((tl >= 7)?3:2);
       b57 = s[ti];
       if (!gabp[*se][b57]) energy -= 2.0;
       else
        { k = (var > 3)?2:((var > 1)?1:0);
          if (bp[cend[k]][b57]) energy += 1.0; }}

  /* R motif at end of tstem */

    if (!RI[s[-1]]) energy -= 2.0;

  /* large tloop */

    if (tl > 9) 
     if (tbondtype > 0x200) energy -= 2.0;

  /* dloop,var,tloop T repeat motif */
  /* present in some nematode D-loop replacement tRNA-Ser genes */

    if (dloop >= 4)
     { k = 1;
       se = aend1;
       while (se < cpos) if (*se++ == Thymine) k++;
       if (k >= dloop)
        { if (var >= 3)
           { se = cend;
             while (se < tpos) if (*se++ == Thymine) k++;
             if (k >= (var + dloop))
              { energy += 3.0;
                se = s + ((tl > 5)?5:tl);
                while (s < se) if (*s++ != Thymine) break;
                if (s >= se) energy += 5.5; }}}}

  /* astem stability  */

    if (ea < -6.1)
     if (tl > 4)
      { if (*s == Thymine)
         if (s[-1] == Guanine)
          if (s[1] == Thymine)
           goto NASI;
        if (ea > -8.3)
         if (*clooppos == Cytosine)
          if (clooppos[1] == Thymine)
           if (cloopend[-2] == Adenine)
            if (cloopend[-1] == Adenine)
             goto NASI;
        energy -= 3.0; }
    NASI:

  /* cstem stability (GC bond count) */

    bondtype = acbondtype + tbondtype;
    if ((cbondtype & 0xf) < 1)
     if ((bondtype & 0xf) < 3) energy -= 1.0;

  /* cloop CTnnnAA motif */

    if (bondtype >= 0x400)
     energy += CX[*clooppos] + TX[clooppos[1]] +
               AXX[cloopend[-1]] + AXX37[cloopend[-2]];
    else
     energy += CX[*clooppos] + TX[clooppos[1]] +
               AX[cloopend[-1]] + AX37[cloopend[-2]];

  /* large dloop */

    if (dloop >= 9)
     { k = tloop - dloop - 4;
       if (k < 0)
        if (bondtype >= 0x1000) energy += (double)k;
         if (dloop >= 12)
          { if (dloop >= 14) energy -= 2.0;
            else
             if (tstem < 6) energy -= ((dloop >= 13)?2.0:1.0); }}

  /* small dloop, small tarm */

    if (dloop <= 10)
     if (tstem < 3) 
      if (ea > -2.6)
        if (tl <= 7)
          if (cgcat >= 4)
            if (gc[*tpos][apos2[-1]])
              if (gc[tpos[1]][apos2[-2]])
               if (gcv > 1.2)
                if ((abondtype & 0xf) > 0)
                  if ((cbondtype & 0xf) > 0)
		   energy += (4.5 + (mtBONDSTAB - 0.5)*(double)(5 - tstem));

  /* global stem damage level */

           i = (int)((bondtype >> 16) & 0xf);
           j = (int)((bondtype >> 12) & 0xf) + i;
           k = (int)((bondtype >> 8) & 0xf);
	   if (tstem > 3)
	    { if ((k > 0) || (tl > 9))
	       if ((j > 0) || (k > 5))
                { n = j + k;
                  if ((s[-1] != Guanine) || (*s != Thymine) || 
                    (s[1] != Thymine) || (tstem < 5))
                   if (n > 4) energy -= 2.0*(double)(n - 4); }}
	   else
	    { n = j + k;
              if (n > 3) energy -= 2.0*(double)(n - 3); }

  /* long tstem with tloop GTT motif */

   if (s[-1] == Guanine)
    if (*s == Thymine)
     if (s[1] == Thymine)
      if (tstem >= 6)
       if (tbondtype < 0x100)
        energy += 1.5;

  /* find exceptions */

           if (energy < tthresh)
            { if (!mtxdetect) goto DN;
              if (incds) goto DN;
              if (energy < (thresh - 13.5)) goto DN;
              if (energy < (tthresh - 13.5)) goto DN;
              if (k > 1)
	       { if (i > 2) goto DN;
                 if (k > 4)
	          if (i > 1) goto DN; }
	      if (nbase > 70) goto DN;
	      if (var > 4)
	       { if (var > 5) goto DN;
	         if (var > tl) goto DN; }
	      if (tstem < 4)
               if ((agcat + cgcat + 2) < (astem + cstem)) 
		goto DN;
              if (tl > 9) goto DN; 
              if (dloop > 13) goto DN;
	      if (!YI[*clooppos]) goto DN;
	      if ((abondtype & 0xf) < 2)
	       { if ((abondtype & 0xf) < 1) goto DN;
                 if (cbondtype > 0x200)
                  if (tbondtype > 0x100)
                   if (abondtype > 0x200)
                    goto DN; }
	      if ((tbondtype & 0xf) < 1)
	       { if ((acbondtype & 0xf) < 1) goto DN;
                 if (acbondtype > 0x200) goto DN; }
              if ((dloop + 19) < tloop) goto DN;
              if (gcv > 5.5) goto DN;
              tgcat = ((tbondtype >> 4) + tbondtype) & 0xf;
	      if ((tgcat + 2) < tstem) goto DN;
	      if (cloop != 7) goto DN;
              if (bp[*cpos][cend[-1]])
                  if (bp[cpos[-1]][*cend])
                   if (bp[cpos[-2]][cend[1]])
		    energy += 2.0;
              if (bondtype < 0x20000)
                  if (thit[nti].stem_energy > -4.6)
                   if (tstem >= 4)
 	            if ((tstem >= 5) || (s[-1] == Guanine))
                     if (stackbp[cpos[1]][cend[-2]])
                      if (stackbp[*cpos][cend[-1]])
                       if (stackbp[cpos[-1]][*cend])
                        { energy += 1.5;
	                  if (s[-1] == Guanine)
                           if (*s == Thymine)
                            if (s[1] == Thymine)
                             energy += 1.0; 
                          if (agcat >= 6) energy += 0.5; }
               if (tc > 0)
                 if (tstem >= 5)
                   if (var < 6)
                    if (var > 2)
		     { if (acbondtype < 0x100) energy += 5.0;
                       else
			if ((abondtype + tbondtype) < 0x100)
			 energy += 3.0;
		        else
                         if (cloopend[-2] == Thymine)
                         if (cloopend[-1] == Thymine)
                         if (dloop > 7)
			 if (tbondtype < 0x100)
                          if (!tt[*tpos][apos2[-1]])
                          if ((agcat+cgcat) >= 10)
                           energy += 13.5; }
	         if (s[-1] == Guanine)
                  if (*s == Thymine)
                   if (s[1] == Thymine)
		    if ((tstem >= 5) || (s[2] == Cytosine))
		     { energy += 1.5;
		       if (tstem >= 5)
                        if (tbondtype < 0x1000)
		         if (s[2] == Cytosine)
                          { if (abondtype < 0x10000)
                             { if (*clooppos == Cytosine)
                                if (clooppos[1] == Thymine)
                                 if (cloopend[-2] == Adenine)
                                  if (cloopend[-1] == Adenine)
                                   energy += 3.0;
			       if (tbondtype < 0x200)
			        if (bondtype < 0x10000)
			         if (tl == 7)
			          if (s[4] == Adenine)
			           energy += 4.0; }}
			 else
			  if (tbondtype < 0x200)
                           if ((tbondtype & 0xf) >= 2)
                             if (*clooppos == Cytosine)
                              if (clooppos[1] == Thymine)
                                 if (cloopend[-2] == Adenine)
                                  if (cloopend[-1] == Adenine)
                                 energy += 1.0; }
		 if (tstem >= 4)
                  if (tbondtype < 0x100)
	           if (cbondtype < 0x200)
		    if (agcat >= 5) energy += 1.5;
	         if (energy > tthresh) energy = tthresh; 
              if (ea > -1.8) energy += 3.0;
              else
	       if (abondtype < 0x60) energy += 1.5;
	       else
		if (acbondtype < 0x200) energy += 0.75;
              if (*clooppos == Cytosine)
                                 if (cloopend[-2] == Adenine)
                                  if (cloopend[-1] == Adenine)
		 { if (tstem >= 5)
                    if (tbondtype < 0x100)
                     if (clooppos[1] == Thymine)
                      { energy += 3.0;
                        if (tstem >= 6) energy += 1.0; }
                     else
                      if (clooppos[1] == Cytosine)
                       energy += 1.0;
		   if (tc >= 2)
		    if (clooppos[1] == Thymine)
                     if (bondtype < 0x1000)
                      if (tstem >= 4)
                        if (var < 6)
                         if (var > 2)
		          energy += 3.0; }
              if (cbondtype < 0x100)
               if (agcat >= 5)
                if (tc > 0)
                 if (clooppos[1] == Thymine)
                  if (YI[*clooppos])
                   if (RI[cloopend[-2]])
                    if (RI[cloopend[-1]])
                     if (tbondtype < 0x100)
                      energy += 4.0;
	             else
                      if (agcat >= 6)
                       if ((tgcat + 1) >= tstem)
                        if (tstem >= 4)
                         energy += 4.0;
              if (bondtype < 0x1000) 
               { energy += 0.5;
                 if (bondtype < 0x200) energy += 0.75; } 
	      if (energy < tthresh) goto DN;
              energy -= (3.0 + 0.9*(energy - tthresh)); }

  /* mammalian cloop motif constraint */

	   if (ds == MAMMAL_MT)
	    { s1 = clooppos;
	      s2 = s1 + cloop;
	      r = *s1++;
	      while (s1 < s2) r = (r << 4) + *s1++;
              if (r != clmotif[0])
               if (r != clmotif[1])
                if (r != clmotif[2])
                 energy -= 5.0; }

  /* remember fully formed D-loop replacement mttRNA gene */
  /* if threshold reached */

       if (energy < thresh) goto DN;
	   te.energy = energy;
       thresh = energy;
       te.ps = apos1;
       te.spacer1 = 0;
       te.spacer2 = 0;
       te.dstem = 0;
       te.dloop = dloop;
       te.cstem = cstem;
       te.cloop = cloop;
       te.anticodon = astem + dloop + cstem + 2;
       te.nintron = 0;
       te.intron = 0;
       te.var = var;
       te.varbp = 0;
       te.tstem = tstem;
       te.tloop = tl;
       te.nbase = astem + dloop + carm + var +
                  2*tstem + tl;
	   tastem = astem;
	   tastem8 = astem8;
	   tastem8d = astem8d;

  /* build fully formed cloverleaf mttRNA genes */

       DN:
       if (dloop < 10) continue;

  /* choose tarm */

          nt = -1;
          nti = -1;
          et = -INACTIVE;
          while (++nt < nth)
	   { tend = thit[nt].end;
             if (tend != apos2) continue;
             e = thit[nt].energy;
	     tpos = thit[nt].pos;
             tstem = thit[nt].stem;

  /* GT motif on tloop */

	     s = tpos + tstem;
	     if (*s == Thymine)
	      if (s[-1] == Guanine)
	       if (tstem >= 5)
	        if (!stackbp[*tpos][tend[-1]])
                 { e += 0.5;
	           if (!bp[tpos[1]][tend[-2]]) e += 0.5; }

  /* large var loop */

             var = (int)(tpos - cend);
             if (var > 5)
              { ev = (double)(var - 5);
		if (tstem < 5) e -= 3.0*ev;
		else e -= (0.5 + 0.5*ev);

  /* allow large var loop if tarm looks nuclear */
  /* (GTTC motif, very large var loop base-pairing) */

                if (var > 9)
                 { if ((thit[nt].bondtype & 0xf) < 1) e -= 1.0;
                   e -= (0.25*(double)(var - 8));
		   if (*s == Thymine)
		    if (s[-1] == Guanine)
		     if (s[1] == Thymine)
		      if (s[2] == Cytosine)
		       e += 4.0;
         if (var > 17)
          { if (var > 25) continue;
            e += 0.5*vloop_stability(cend,var,&varbp); }}}

  /* small var loop */

             if (var < 3)
              { if (tstem > 5)
	         if (s[-1] != Guanine)
		  e -= 0.5;
                if (var < 2)
	         { if (var < 1)
                    { if (var < 0) continue;
                      if (tstem < 4)
                       if (thit[nt].stem_energy < -4.0)
                        continue; }
                   e -= 3.0; }}
             if (e > et)
              { et = e;
                nti = nt; }}

          if (nti < 0) continue;
          tpos = thit[nti].pos;
          tstem = thit[nti].stem;
          tl = thit[nti].loop;
          tarm = 2*tstem + tl;
          var = (int)(tpos - cend);
          b48 = tpos[-1];
          tbondtype = thit[nti].bondtype;
	  bondtype = acbondtype + tbondtype;
	  ti = (int)(((bondtype >> 16) & 0xf) + ((bondtype >> 12) & 0xf) +
               ((bondtype >> 8) & 0xf));

  /* choose darm */

          nd = -1;
          ndi = -1;
          ed = -INACTIVE;
          while (++nd < ndh)
           { dl = dhit[nd].loop;
             dstem = dhit[nd].stem;
             darm = 2*dstem + dl;
             dpos = dhit[nd].pos;
             e = dhit[nd].energy;

  /* spacing between astem,darm,carm */

             spacer1 = (int)(dpos - aend1);
             spacer2 = (int)(cpos - dpos) - darm;
             if (spacer1 < 2)
              { if (spacer1 < 1) continue;
		if (dstem < 3) continue;
                if (dl > 12) e -= 2.0;
		if (astem < 7) e -= 1.0;
                if (spacer2 != 2)
                  { if (spacer2 < 1) continue;
                    if (spacer2 > 2) continue;
		    if ((abondtype & 0xf) < 1)
		     if ((dhit[nd].bondtype & 0xf) < 1)
		      e -= 0.5;
		    if (var > 7) e -= 1.0;
		    if (dl > 12) e -= 1.0;
		    if (cloop != 7) e-= 2.0;
                    if (cstem < 6) e -= 3.6;
                    else e -= 0.5; }
                 else
                  { if (cstem > 5) continue;
                    s = cpos;
                    se = cend-1;
                    while (!bp[*s][*se])
                     { s++;
                       se--; }
                    if (!stemterm[s[-1]][se[1]]) e -= 0.5;
                    e -= 0.8; }}
              else
               { if (spacer1 > 2)
		  { if (spacer1 > 3) continue;
		    if (dstem > 4) continue;
		    if (dstem < 3) continue;
		    if (tl > 15) continue;
		    if (astem < 7) e -= 1.0;
		    if (ti > 4) e -= 1.0;
		    if (cloop != 7) e-= 2.0;
		    if (tbondtype > 0x2000)
		     if (!RI[tpos[tstem-1]]) e -= 2.0;
	            e -= 1.0;
		    if (spacer2 != 1) e -= 0.5;
		    else
		     if (dhit[nd].bondtype < 0x100)
		      if (var >= 3)
                       if (var <= 5)
                        if (tstem >= 3)
			 { e += 1.0;
                           if (agcat >= 5)
                            if (wcbp[*aend1][*apos2])
                             if (!bp[aend1[-1]][*apos2])
		              if (bp[b48][dpos[dstem+1]])
                               e += 0.5; }
		    }
                 if (spacer2 > 1)
                  { if (spacer2 > 2) continue;
                    if (astem < 7)
                     if (spacer1 == 2)
                      e -= 1.0;
		    if (cloop != 7) e -= 2.0;
		    if (ea < -5.8) e -= 2.0; 
                    e -= 2.5; 
		    if (bp[b48][dpos[dstem+1]]) 
                     { if (dhit[nd].bondtype < 0x1000)
                        if (wcbp[dpos[1]][dpos[darm-2]])
                         if (wcbp[dpos[2]][dpos[darm-3]])
                          if (var < 6)
                           if (dl > 3)
                            e += 2.0; }
                    else e -= 1.0;
		     }
                 else
		  if (spacer2 < 1)
		   { if (spacer2 < 0) continue;
                     if (var > 6) continue;
		     if (dstem > 4) continue;
		     if (dhit[nd].stem_energy < -4.3) continue;
                     if (astem < 7)
                      if (spacer1 == 2)
                       e -= 1.0;
		     if (cloop != 7) e-= 2.0;
		     e -= mtBONDSTAB; }
                 if (cstem > 5)
                  if ((!gt[*cpos][cend[-1]]) || astem8) e-= mtBONDSTAB;  }

  /* very large or very small dloop */

        if (dl < 3) e -= 2.0;
        if (dl > 11)
         { if (dl > 14) e -= 2.0;
	       else
	        if (dl > 13)
	         { if (dhit[nd].bondtype >= 0x100) e -= 2.0;
	           else e -= 1.0; }
	        else
	         if (dl > 12)
	          { if (dhit[nd].bondtype >= 0x1000) e -= 2.0;
	            else e -= 1.0; }
             else
	           if (dhit[nd].bondtype >= 0x10000) e -= 2.0; }

  /* tertiary interactions in darm */

             b8 = dpos[-2];
             b9 = dpos[-1];
             if (dl > 2)
               { if (dl > 5)
                  if (!stackbp[dpos[dstem+1]][b48]) e -= 1.0;
                 if (!stackbp[b8][dpos[dstem]]) e-= 0.25;
                 if (!stackbp[b8][dpos[dstem+dl-1]]) e -= 0.25; }
              if (!bp[b9][dpos[2]])
               if (!bp[b9][dpos[darm-3]])
                e -= 1.0;

  /* TR motif at b8-9 */

              if (RI[b9])
               { if (b8 == Thymine)
                  if (spacer1 == 2)
                   if (ti < 6)
                    if (((bondtype & 0xf) > 2) || (bondtype < 0x1000) ||
	                ((tbondtype < 0x100) && (tstem > 3)))
                     if ((cbondtype & 0xf) < 5)
                    if (stembp[dpos[1]][dpos[darm-2]])
                    if (var < 6)
                     if (var > 2) e += 1.5;
		     else
                      if (tstem > 3)
			   if (cloopend[-2] == Adenine)
		        e += 1.5; }
              else
               { e -= 1.0;
                 if (b9 == Thymine)
                  if (spacer1 == 2) e -= 2.0; }
              if (e > ed)
               { ed = e;
                 ndi = nd; }}

          if (ndi < 0) continue;
          energy = 100.0 + ec + ed + et;
          dl = dhit[ndi].loop;
          dstem = dhit[ndi].stem;
          darm = 2*dstem + dl;
          dpos = dhit[ndi].pos;
          dbondtype = dhit[ndi].bondtype;
          spacer1 = (int)(dpos - aend1);
          spacer2 = (int)(cpos - dpos) - darm;
          b8 = dpos[-2];

  /* tertiary structure interaction between tloop and dloop */

          if (tl >= 3)
           if (dl >= 4)
            { di = (dl < 7)?(darm-dstem-2):(darm-dstem-3);
              ti = (tl < 9)?(tstem+2):((tl < 13)?(tstem+3):(tstem+5));
              if (ggbp[dpos[di]][tpos[ti]])
               if (ggbp[dpos[di-1]][tpos[ti-1]])
                { energy += 2.0;
                  if (spacer1 != 2)
                   if (spacer2 != 2)
                    if (dstem < 4)
                     if (tl > 7)
                      if (bp[dpos[di+1]][tpos[ti+1]])
                       energy += 4.0;
                  if (ea > -2.5)
                   if (wcbp[dpos[1]][dpos[darm-2]])
                    if (wcbp[dpos[2]][dpos[darm-3]])
                     energy += 3.0; }
              if (tl > 10)
               if (dl > 10)
                energy -= 1.0; }
           else
            if (dl == 3)
             if (wcbp[dpos[dstem+1]][b48]) energy += 1.0;

  /* small darm and tarm */

         if (tloop <= 18)
          if (tarm <= 13)
           if (dl <= 8)
            if (spacer1 == 2)
             if (spacer2 == 1)
              if (abondtype < 0x1000)
               if (tbondtype < 0x100)
                if (dbondtype < 0x200)
                { et = (mtBONDSTAB - 0.5)*(double)(5 - tstem) + 
			0.1*(double)(7-tl);
                  ed = mtBONDSTAB*(double)(4 - dstem);
                  energy += (0.8*(et + ed)); }

  /* GTTC motif on tloop */

          s = tpos + tstem;
          if (tl < 5)
           if (tl < 2) energy += G[s[-1]];
           else
           { et = (G[s[-1]] + T[*s] + T[s[1]]);
             if (tl > 3)
             if (bp[*s][s[tl-1]])
              { e = (G[*s] + T[s[1]] + T[s[2]]);
                if (e > et) et = e; }
             if (tstem < 5)
              { e = (G[s[-2]] + T[s[-1]] + T[*s] + C[s[1]]);
                if (e > et) et = e; }
             energy += et; }
         else energy += (G[s[-1]] + T[*s] + T[s[1]] + C[s[2]]);

  /* long astem */

         if (astem8)
          if (bp[apos1[0]][apos2[6]])
           if (bp[apos1[1]][apos2[5]])
            if (bp[apos1[2]][apos2[4]])
             if (bp[apos1[3]][apos2[3]])
              energy += hbem[apos1[-1]][apos2[7]];

  /* false positive supression */

         if (!RI[cend[0]]) energy -= 1.0;
         if (!RI[cpos[-1]]) energy -= 1.0;
	 if (tarm < (var + 3)) energy -= 2.0;
         if (gcv < 1.5)
          if (dbondtype > 0x10000) energy -= 2.0;
	 if (tarm > 27)
          { energy -= 1.0;
            if (spacer2 != 1) energy -= 1.0; }
	 if (dstem < 3)
	  { if (var > 5) energy -= 1.0;
            if (tloop > (dloop + 8)) energy -= 0.5; }
	 if (b8 != Thymine)
          if (dl > 3)
           if (dbondtype > 0x100)
	    if ((b8 == Cytosine) ||
               (dbondtype > 0x10000))
             if (*clooppos != Cytosine)
              if (!wcbp[dpos[dstem+1]][b48])
               energy -= 1.0;


  /* high GC false positive suppression */

	 if (gcv >= 5.1)
	  { if ((abondtype & 0xf) >= 4)
             { s1 = apos1;
               s2 = apos2 + astem;
               n = 0;
               while (--s2 >= apos2)
                if (gc[*s1++][*s2])
                 { if (++n >= 4)
		    { energy -= 2.0;
                      break; }}
                else n = 0; }
	    if ((dbondtype & 0xf) >= 4) energy -= 3.0;
	    if ((cbondtype & 0xf) >= 5) energy -= 3.5;
	    if ((tbondtype & 0xf) >= tstem) energy -= 4.0; }

  /* global stem damage level */

	 tc = tstem + dstem;
	 dtbondtype = dbondtype + tbondtype;
         mabondtype = dtbondtype + cbondtype;
         bondtype = acbondtype + dtbondtype;
         if (bondtype < 0x100) energy += 0.5;
	 if ((dtbondtype & 0xf) < 1)
          { energy -= 1.0;
            if (tc >= 10) energy -= 2.0;
	    if ((bondtype & 0xf) < 3)
             if (nbase > 75) energy -= 1.0; }
         i = (int)((bondtype >> 16) & 0xf);
         j = (int)((bondtype >> 12) & 0xf) + i;
         k = (int)((bondtype >> 8) & 0xf) + j;
	 ti = (tc > 6)?5:((tc > 5)?4:3);
         if (k > ti)
          { ev = (double)(k - ti);
            energy -= 0.5*ev;
            if (cbondtype > 0x10000)
             if (tstem < 5)
              energy -= ev;
            if (i > 0)
             if (k > 8)
              energy -= 1.5*(double)(k - 8); }

  /* low GC false positive supression */

         if (gcv < 3.5)
          if ((bondtype & 0xf) < 2)
           { if ((bondtype & 0xf) < 1) energy -= 1.0;
              if (dl > 3)
               if (var > 2)
                if (!wcbp[dpos[dstem+1]][b48]) energy -= 1.0; }

  /* small variable loop */

         if (var < 3)
	  { if (dloop > 18)
	     { if (dloop > (tloop + 2)) energy -= 1.0;
	       if (tloop > 20)
                if ((((dtbondtype >> 4) + dtbondtype) & 0xf) < 6)
	         energy -= 2.0; }
	    if (astem < 7) 
             { energy -= 1.0;
               if (agcat >= 5)
	        if (bondtype < 0x300)
	         if (gcv > 1.2)
                  if (gcv < 5.0)
                   energy += 2.0; }}
	 else

  /* NNNNNAA cloop */

          if (cloopend[-2] == Adenine)
           if (cloopend[-1] == Adenine)
            if (spacer1 > 1)
            if ((dbondtype < 0x2000) || (dloop > mt_DRLmaxlength))
             { if (abondtype < 0x100) energy += 1.0;
               else
		if (cbondtype < 0x100) energy += 1.0;
	        else
                 if (tstem >= 5)
		  if (tbondtype < 0x100)
                   { energy += 1.0;
                     if (*clooppos == Cytosine)
		      if (clooppos[1] == Thymine)
                       if (dbondtype < 0x100)
		       energy += 0.5;

                     if (cgcat >= 3)
                     if ((tbondtype & 0xf) > 0)
                     if (ggbp[dpos[dstem+1]][b48])
                     if (wcbp[dpos[1]][dpos[darm-2]])
                     if (tl < 10)
                     if (spacer1 == 2)
	              if (spacer2 == 1)
                       if (dl > 2)
                        if (var >= 2)
                         if (var < 6)
                          { if (agcat >= 6) energy += 3.0;
			    else
			     if (agcat >= 5)
                              if (cgcat >= 4)
                               if (dbondtype < 0x100)
                                if (*s == Thymine)
                                 if (s[-1] == Guanine)
                                  if (s[1] == Thymine)
                                   energy += 3.0; }}}

  /* large tloop */

         if (tl > 12)
	  { if (tbondtype > 0x10000) energy -= 2.0;
            if (agcat < 5)
             if (spacer1 != 2)
              if (spacer2 != 1)
               energy -= 1.0; }

  /* find exceptions */

	 if (energy < dtthresh)
      { if (!mtxdetect) continue;
        if (incds) continue;
        if (energy < (thresh - 12.0)) continue;
        if (energy < (dtthresh - 12.0)) continue;
	    if (nbase > 75) continue;
	    if (dstem > 4) continue;
	    if (dstem < 3) continue;
	    if (astem < 7)
	     if (acbondtype > 0x21000)
              continue;
            if (var > 5)
             { if (var > 6) continue;
               if (tarm < 12) continue; }
	    if (gcv <= 1.2)
             { if (gcv < 0.9) continue;
               if ((mabondtype & 0xf) < 1) continue; }
	    if (tl > 9)
             { if (tl > 13) continue;
               if (!wcbp[dpos[1]][dpos[darm-2]])
                continue; }
	    if (dl > 7)
         { if (bondtype > 0x20000)
           if (dloop > (tloop + 4)) continue;
	       if (dl > 10)
            { if (dl > 12)
               if (abondtype > 0x1000)
		        continue;
	          if (tbondtype > 0x200) continue;
              if (tt[*tpos][apos2[-1]]) continue;
              if (var > 5) continue;
	          if (dloop > (tloop + 8))
		       if (bondtype > 0x10000) continue;
		      if (astem < 7) continue; }}
	    if (RI[clooppos[1]]) continue; 
	    b9 = dpos[-1];
	    if (cstem >= 6)
	     { if (cbondtype > 0x200) continue;
               if (var < 3) continue; 
	       if (YI[b9]) continue; }
	    if (cloop != 7) continue;
            if (ds == MAMMAL_MT) continue;
            if (mabondtype < 0x400)
             { if ((b8 == Thymine) || (mabondtype < 0x300))
		if (ea < -5.45)
	         if (chit[nc].stem_energy > -3.2)
		  if (dbondtype < 0x200)
                  if (spacer1 > 1)
                   if ((spacer2 == 1) || (mabondtype < 0x100))
                    if ((spacer1 < 3) || (tstem > 3) ||
			(tbondtype < 0x100))
                     if ((spacer1 < 3) ||
		      ((var > 2) && (var < 6) && (tbondtype < 0x2000) &&
		       (tl < 10)))
                    if (dstem < 5)
                   if (var >= 2)
 	            if (dl > 2)
                     if (tl < 15)
                      if ((b8 != Cytosine) ||
		        (*clooppos == Cytosine))
                     if (RI[b9])
                      if (*clooppos != Adenine)
                       if (clooppos[1] == Thymine)
			if (RI[cloopend[-2]])
                         { s1 = apos1;
                           s2 = apos2 + astem;
                           n = 0;
                           while (--s2 >= apos2)
                            if (wcbp[*s1++][*s2])
                             { if (++n >= 3) break; }
                            else n = 0;
	                   if (n >= 3)
			    { energy += 3.0;
                              if ((abondtype & 0xf) > 0) energy += 2.0;
                              if (bp[dpos[dstem+1]][b48])
                                if (wcbp[dpos[1]][dpos[darm-2]])
                                 if (var <= 5)
                                  energy += 1.0; }
			   if (dtbondtype < 0x200)
                           if (agcat < 2)
                            if (wcbp[dpos[dstem+1]][b48])
                             if (wcbp[dpos[1]][dpos[darm-2]])
                              if (wcbp[dpos[2]][dpos[darm-3]])
			      if (gcv > 1.2)
                              if (var <= 5)
			       if (tstem >= 3)
			        if (dstem >= 3)
                                 if (tl > 3)
                                  if (tl < 9)
                                   if (dl < 9)
                                    if (spacer1 == 2)
				     energy += 10.0; }
               if ((tbondtype & 0xf) > 0)
               if (mabondtype < 0x300)
	        { if (mabondtype < 0x100)
                   { if ((spacer1 < 3) || (tstem > 2))
                    if (var > 0)
                     if (YI[*clooppos])
                      if ((spacer2 > 0) || (clooppos[1] == Thymine))
	               energy += 2.5; }
	        else
                 if ((dbondtype & 0xf) > 0)
	          if (b9 != Cytosine)
                   if (var <= 7)
		    if (spacer2 == 1)
                     if (tarm < 22)
                      if (gcv > 1.2)
                     if (dstem >= 4)
		      { if (tstem >= 5)
                       energy += 5.0;
	              else
	               if (tstem >= 3)
	                if (tbondtype < 0x100)
	                 energy += 1.0; }
	             else
		      if (tstem >= 5)
                        energy += 1.0; }
             else
              if ((dbondtype & 0xf) > 0)
              { if (tstem >= 5)
		 if (s[-1] == Guanine)
		  if (*s == Thymine)
		   if (s[1] == Thymine)
                    if (*clooppos == Cytosine)
                     if (clooppos[1] == Thymine)
			  if (cloopend[-2] == Adenine)
                       if (cloopend[-1] == Adenine)

			energy += 1.0;
	         if (bondtype < 0x1000)
	          if (cbondtype < 0x100)
	           energy += 1.0; }}
			
	    if (tstem >= 5)
		 if (*clooppos == Cytosine)
		 {
             if (dl > 3)
             if (dtbondtype < 0x200)
              if ((tbondtype & 0xf) > 0)
              if (clooppos[1] == Thymine)
               {
	        if (clooppos[2] == Thymine)
		  if (clooppos[3] == Adenine)
              if (clooppos[4] == Cytosine)
              if (clooppos[5] == Adenine)
			  if (cloop == 7)
		      energy += 0.5;

                if (cgcat >= 4)
		 if (wcbp[dpos[1]][dpos[darm-2]])
                  if (bp[dpos[dstem+1]][b48])
                     if (tl < 10)
                      if (var < 6)
	               if (spacer1 == 2)
	                if (spacer2 == 1)
		        if (dstem >= 3)
		         energy += 3.0;

	       }
			
			
            if (clooppos[1] == Cytosine)
              if (clooppos[2] == Cytosine)
		  if (clooppos[3] == Adenine)
	        if (clooppos[4] == Thymine)
  		  if (s[-1] == Guanine)
		   if (*s == Thymine)
		    if (s[1] == Thymine)
		     energy += 1.0;
		}	

		
		if (RI[b9])
		{	
		if (b8 == Thymine)
		{
		if (clooppos[1] == Thymine)
                {
		if (cloopend[-2] == Adenine)
		{
		
		if (wcbp[dpos[1]][dpos[darm-2]])
		{
                if (*clooppos == Cytosine)
	        {
	         if (abondtype < 0x200)
		   energy += 1.0;

                 if (bondtype < 0x10000)
                  if (dtbondtype < 0x200)
                   if (agcat >= 3)
                    if (cgcat >= 4)
                     if (tl < 10)
                      if (var < 6)
	               if (spacer1 == 2)
	                if (spacer2 == 1)
                          if (tstem >= 3)
		            energy += 3.0;
		}

			
                if (tstem >= 5)
		 if (s[-1] == Guanine)
		  if (*s == Thymine)
		   if (s[1] == Thymine)
                    { energy += 1.0;
                      if (tl >= 5)
	               if (spacer1 == 2)
	                if (spacer2 == 1)
                         if (tbondtype < 0x100)
                          if (wcbp[dpos[dstem+1]][b48])
                           energy += 3.0; }

                if (tstem >= 3)
                if (tl < 10)
	        if (spacer1 == 2)
	        if (spacer2 == 1)
		 if (RI[cloopend[-1]])
                  if (dl > 2)
                   if (var >= 2)
                  if (var < 6)
                  if (ggbp[dpos[dstem+1]][b48])
		  {
                   if (dtbondtype < 0x100)
                    { energy += 3.5;
                      if ((bondtype & 0xf00) == 0)
                       if (*clooppos == Cytosine)
                         energy += 1.5; }


                   if (bondtype < 0x10000)
	            if (tstem > 2)
                     if (tbondtype < 0x200) energy += 2.5;

		   if (abondtype < 0x100)
		    if (wcbp[dpos[2]][dpos[darm-3]])
                     energy += 3.0;

                   if (tbondtype < 0x100)
                   if (agcat >= 6)
                   if (tstem >= 5)
                   if ((tbondtype & 0xf) > 0)
                   if (RI[cloopend[-1]])
                   if (cgcat >= 4)
                    energy += 2.0;
		  }
		  else
                   if (!ggstembp[*tpos][apos2[-1]])
                    if (wcbp[dpos[dstem+1]][*tpos])
                     energy += 1.5;
        }
			
	    if ((abondtype & 0xf) < 1)
	    if (abondtype < 0x100)
             if (gcv > 1.2)
	      if (dl > 3)
                if (bp[dpos[dstem+1]][b48])
                if (spacer1 == 2)
                 if (spacer2 == 1)
             if (*clooppos == Cytosine)
		      energy += 5.0;

	    if (cbondtype < 0x100)
             if (tbondtype < 0x100)
              if (tstem >= 3)
	      if (dl > 3)
               if (var < 6)
                if (bp[dpos[dstem+1]][b48])
                if (spacer1 == 2)
                 if (spacer2 == 1)
		  energy += 2.5;

        }

        if (stembp[dpos[dstem+1]][b48])
		{

             if (*clooppos == Thymine)
                if (cloopend[-2] == Guanine)
              if (clooppos[2] == Guanine)
	       if (clooppos[3] == Thymine)
	        if (clooppos[4] == Guanine)
		  if (dl > 2)
              energy += 1.0;
				
	    if (cbondtype < 0x100)
             if (dbondtype < 0x10000)
                  if (wcbp[dpos[1]][dpos[darm-2]])
		    if (var < 6)
		     if (tstem >= 3)
			if (gcv >= 1.2)
			  if (dl > 3)
                energy += 1.0;

	    if (tstem >= 5)
             if (dtbondtype < 0x200)
              if (*clooppos == Cytosine)
                if (spacer1 == 2)
                 if (spacer2 == 1)
		 if (RI[cloopend[-2]])
                  energy += 0.5;
		}	

		}		
				

	     if (tstem > 2)
              if (tarm < 28)
                if (spacer1 == 2)
                 if (spacer2 == 1)
                  if (dl > 3)
                   if (j < 1)
		    if (k > ti)
                   if (ggstembp[dpos[dstem+1]][b48])
			energy += 2.5;
			
            if (dtbondtype < 0x100)
              if ((tbondtype & 0xf) > 0)
               if (bp[dpos[dstem+1]][b48])
                 if (b9 == Adenine)
                  if ((dbondtype & 0xf) > 0) energy += 2.0;
	          else
                   if (spacer2 == 1)
	            energy += 0.5;
	    if (cloopend[-2] == Adenine)
	     if (cloopend[-1] == Adenine)
              if (cbondtype < 0x2000)
               if (spacer1 > 1)
 	       if (dl > 2)
                if (var < 6)
                  energy += 0.75;
				
		
		}
		
            if (var > 2)
             if (dl > 2)
              {
              if (cbondtype < 0x200)
               if (((mabondtype & 0xf) > 3) || (bondtype < 0x1000))
		{
                if (bp[dpos[dstem+1]][b48])
                 energy += 1.0;

                if (cbondtype < 0x100)
		if (dbondtype < 0x100)
	        if (bp[b8][dpos[dstem]])
                if (bp[b8][dpos[darm-dstem-1]])
                if (var < 6)
                if (tstem >= 3)
                if (tl < 10)
	        if (spacer1 == 2)
	        if (spacer2 == 1)
                if (clooppos[1] == Thymine)
		 if (cloopend[-2] == Adenine)
                  energy += 3.0;
		}

              if (clooppos[1] == Thymine)
		   if (RI[cloopend[-2]])
		{
                if (*clooppos == Cytosine)
		{

	      if (dtbondtype < 0x200)
                if (agcat >= 3)
                if (cgcat >= 4)
                if (var < 6)
                if (tstem >= 3)
                if (tl < 10)
	        if (spacer1 == 2)
	        if (spacer2 == 1)
                 { if (abondtype > 0x20000)
                    if (bp[dpos[dstem+1]][b48])
		     energy += 7.0;
                   if (agcat >= 6) energy += 2.0; }
	
	        if ((bondtype & 0xf00) == 0)
		 if (gcv > 5.0)
  		  if (s[-1] == Guanine)
		   if (*s == Thymine)
                    if (tstem >= 5)
                    if (var < 6)
                    if (tl < 10)
	            if (spacer1 == 2)
	            if (spacer2 == 1)
		     energy += 2.0;

		if (abondtype < 0x100)
                 if (cbondtype < 0x10000)
                  if (bp[dpos[dstem+1]][b48])
                   if (cgcat >= 4)
                    if (tstem >= 3)
                     if (var < 6)
                      if (tstem >= 3)
                       if (tl < 10)
	                if (spacer1 == 2)
	                 if (spacer2 == 1)
			  energy += 1.5;
	        }


	      if (dtbondtype < 0x100)
               if (agcat >= 4)
               if (cgcat >= 4)
                if (var < 6)
                if (tstem >= 3)
                if (tl < 10)
                 { if (spacer1 == 2)
	            { if (abondtype < 0x3000)
                       if (stackbp[dpos[dstem+1]][b48])
		        energy += 3.0;
		      if (b8 == Thymine)
  		       if (s[-1] == Guanine)
		        if (*s == Thymine)
		         if (s[1] == Thymine)
		          energy += 3.5; }
                   if (agcat >= 6)
		    if (YI[*clooppos])
  		     if (s[-1] == Guanine)
		      if (*s == Thymine)
                       if ((dtbondtype & 0xf) > 0)
                        energy += 3.0; }

              if (mabondtype < 0x10000)
               if (dtbondtype < 0x400)
                if (agcat >= 5)
                 if (cgcat >= 3)
                  if (tl < 10)
                   if (var < 6)
	            if (spacer1 == 2)
	             if (spacer2 == 1)
		     {
		
                     if (dtbondtype < 0x200)
                      if (cbondtype < 0x300)
                      if (bondtype < 0x10000)
                       if (tstem >= 3)
		         energy += 1.0;

		     if (tstem >= 5)
                      if (s[-1] == Guanine)
		       energy += 4.0;
		     }
		}
	      }
		
		}
		else
                 if (bondtype < 0x10000)
                 if (mabondtype < 0x500)
                 if (dbondtype < 0x100)
		 if (b8 == Thymine)
                 if (agcat >= 4)
		 if (clooppos[1] == Thymine)
		 if (cloopend[-2] == Adenine)
		 if (cloopend[-1] == Adenine)
	         if (spacer1 == 2)
	         if (spacer2 == 1)
                 if (dstem >= 3)
                 if (tstem >= 5)
                 if (dl > 2)
                 if (tl < 10)
                 if (var < 6)
		  energy += 7.0;

           if (agcat >= 5)
            if (cgcat >= 4)
             if ((acbondtype & 0xf) >= 3)
	     { if (tbondtype < 0x100)
                if ((dbondtype & 0xf) > 0)
                if ((((dbondtype >> 4) + dbondtype) & 0xf) >= 3)
                if (wcbp[dpos[dstem+1]][b48])
                if (b8 == Thymine)
                  if (RI[b9])
		   if (clooppos[1] == Thymine)
                   if (YI[*clooppos])
		   if (RI[cloopend[-2]])
		   if (RI[cloopend[-1]])
	           if (spacer1 == 2)
	           if (spacer2 == 1)
                   if (dl > 2)
                   if (tl < 10)
                   if (var < 6)
                   if (var > 2)
                    energy += 6.0;
               if (cgcat >= 5)
                if (abondtype < 0x10000)
                 if (bp[dpos[dstem+1]][b48])
		   if (clooppos[1] == Thymine)
                   if (YI[*clooppos])
		   if (RI[cloopend[-2]])
                   if (dl > 2)
                   if (tl < 10)
                   if (var < 6)
                   if (var > 2)
                    energy += 6.0; }

	    if (energy >= dtthresh)
             energy -= (0.9*(energy - dtthresh) + 5.0);
	    else continue; }

  /* remember fully formed mttRNA gene if threshold reached */

            if (energy < thresh) continue; 
	        te.energy = energy;
            thresh = energy;
            te.ps = apos1;
            te.spacer1 = spacer1;
            te.dstem = dstem;
            te.dloop = dl;
            te.spacer2 = spacer2;
            te.cstem = cstem;
            te.cloop = cloop;
	        te.var = var;
	        te.varbp = (var > 17)?varbp:0;
            te.tstem = tstem;
            te.tloop = tl;
         k = astem + spacer1 + darm + spacer2;
         te.anticodon = k + cstem + 2;
         te.nintron = 0;
         te.intron = 0;
         te.nbase = k + carm + var + 2*tstem + tl;
	 tastem = astem;
	 tastem8 = astem8;
	 tastem8d = astem8d;
         } }

  /* for highest energy mttRNA gene */ 
  /* decide astem length, look for NCCA acceptor tail */
  /* and calculate total length */

    if (te.ps)
     {    apos2 = te.ps + te.nbase;
          if (extastem)
           if (tastem8d)
            { te.astem1 = 8;
              te.astem2 = 8;
              te.ps--;
	      te.nbase++;
	      te.anticodon++;
              as = aatail(apos2+8,&aext,sw); }
           else
            { te.astem1 = tastem;
              te.astem2 = tastem;
              as = aatail(apos2+tastem,&aext,sw);
              if (tastem8)
               { as8 = aatail(apos2+8,&aext8,sw);
                 if (as8 >= as)
                  { te.ps--;
	            te.nbase++;
	            te.anticodon++;
                    te.astem1 = 8;
                    te.astem2 = 8;
                    as = as8;
                    aext = aext8; }}}
          else
           { te.astem1 = tastem;
             te.astem2 = tastem;
             as = aatail(apos2+tastem,&aext,sw); }
          if (as < 2) aext = 1;
          te.nbase += te.astem2;
          nbasefext = te.nbase + ASTEM2_EXT;
          te.nbase += aext;

  /* store mttRNA gene if there are no */
  /* higher energy overlapping mttRNA genes */

          te.start = (long)(te.ps - seq);
          if (tn = find_slot(d,&te,&nts,sw))
           { base_copy3(te.ps,te.seq,nbasefext);
             base_copy3(te.ps,te.eseq,nbasefext);
             te.aatail = aext;
             *tn = te;  }}
     }
  return(nts); }



int tmopt(data_set *d,
          trna_loop *th, int tarm, double the,
          trna_loop *ahit, int nah,
          int nts,int *seq, csw *sw)
{ int r,na,nr,nrh,ibase,flag,as,aext,nbasefext;
  int *s,*v,*s1,*s2,*sa,*sb,*se,*sf,*ps,*tpos,pseq[MAXETRNALEN+1];
  static int gtem[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
  static double A[6] = { 6.0,0.0,0.0,0.0,0.0,0.0 };
  static double Ar[6] = { 10.0,0.0,0.0,0.0,0.0,0.0 };
  static double Cr[6] = { 0.0,10.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double Ga[6] = { 0.0,0.0,7.0,0.0,0.0,0.0 };
  static double K[6] = { 0.0,0.0,6.0,6.0,0.0,0.0 };
  static double Tr[6] = { 0.0,0.0,0.0,10.0,0.0,0.0 };
  double e,energy,penergy,tenergy,aenergy,athresh,cthresh,cathresh;
  static double bem[6][6] =
   { { -1.072,-0.214,-1.072, ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 3.000,-1.072, 0.000, 0.000 },
     { -1.072, 3.000,-1.072, 1.286, 0.000, 0.000 },
     {  ATBOND,-1.072, 1.286,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static trna_loop rhit[NH];
  gene te,*tn;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,0,0,0,13,8,0,28,0,0,3,0,5,7,
     tmRNA,0.0,0,0,0 };
  tpos = th->pos;
  flag = 0;
  te.energy = sw->tmrnathresh;
  athresh = sw->tmathresh;
  cthresh = sw->tmcthresh;
  cathresh = sw->tmcathresh;
  s = tpos + tarm + 4;
  v = tpos + th->stem - 10;
  energy = K[*v] + G[v[1]] + A[v[2]];
  e = K[v[1]] + G[v[2]] + A[v[3]];
  if (e > energy) energy = e;
  if (energy < 18.0) energy = 0.0;
  tenergy = Tr[*s]+Cr[s[1]]+Cr[s[2]]+Ar[s[3]] + energy + 1.59*the;
  nrh = find_resume_seq(tpos-MAXTPTSDIST,TPWINDOW,rhit,NH,sw);
  nr = -1;
  while (++nr < nrh)
   { ps = rhit[nr].pos;
     penergy = tenergy + rhit[nr].energy - 0.001*((double)(tpos - ps));
     if (rhit[nr].stem < 24) penergy -= 15.0;
     na = -1;
     while (++na < nah)
      { aenergy = ahit[na].energy;
        if (aenergy < athresh) continue;
        t.ps = ahit[na].pos;
 if (t.ps < (ps - MAXTPDIST)) continue;
 if (t.ps > (ps - MINTPDIST)) break;
 energy = -INACTIVE;
 sa = t.ps + t.astem1;
        for (sb=sa+9, se=sb+t.cstem; sb <= (sa+16); sb++,se++)
  for (sf = tpos-3; sf >= (tpos-7); sf--)
   { s1 = sb;
            s2 = sf;
     e = bem[*s1++][*--s2];
     while (s1 < se) e += bem[*s1++][*--s2];
     if (e > energy)
      { energy = e;
        t.var = (int)(tpos - sf);
        t.dloop = (int)(sb - sa); }}
        if (energy < cthresh) continue;
        energy += aenergy;
        if (energy < cathresh) continue;
        sb = sa + 3;
        sf = sa + 7;
        r = gtem[*sb++];
        while (sb < sf)
         { r = (r >> 4) + gtem[*sb++];
           if ((r & 3) == 2)
            { energy += 14.0;
              break; }}
        t.energy = penergy + Ga[t.ps[1]] + Ga[t.ps[2]] + energy;
        if (t.energy > te.energy)
         { flag = 1;
    t.tstem = th->stem;
    t.tloop = th->loop;
           t.tps = (int)(ps - t.ps);
           t.tpe = t.tps + rhit[nr].stem;
    ibase = (int)(tpos - t.ps);
           t.nintron = ibase - t.var - 2*t.cstem -
                       t.dloop - t.astem1;
           t.nbase = ibase + tarm + t.astem2 - t.nintron;
           te = t; }}}
  if (flag)
   { te.start = (long)(te.ps - seq);
     s = te.ps + te.nbase + te.nintron;
     as = aatail(s,&aext,sw);
     nbasefext = te.nbase + ASTEM2_EXT;
     te.nbase += aext;
     tn = find_slot(d,&te,&nts,sw);
     if (tn)
      { te.intron = te.astem1 + te.dloop + te.cstem;
        te.asst = 0;
        base_copy3(te.ps,te.eseq,nbasefext+te.nintron);
        remove_intron(te.ps,pseq,nbasefext,
                      te.intron,te.nintron);
        base_copy3(pseq,te.seq,te.nbase);
       te.aatail = aext;
       *tn = te; }}
  return(nts); }


int tmopt_perm(data_set *d,
          trna_loop *th, int tarm, double the,
          trna_loop *ahit, int nah,
          int nts, int *seq, csw *sw)
{ int r,na,nr,nrh,flag,as,aext;
  int *s,*v,*s1,*s2,*sa,*sb,*se,*sf,*ps,*apos,*tpos;
  static int gtem[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
  double e,energy,penergy,tenergy,aenergy,athresh,cthresh,cathresh;
  static double A[6] = { 6.0,0.0,0.0,0.0,0.0,0.0 };
  static double Ar[6] = { 10.0,0.0,0.0,0.0,0.0,0.0 };
  static double Cr[6] = { 0.0,10.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double Ga[6] = { 0.0,0.0,7.0,0.0,0.0,0.0 };
  static double K[6] = { 0.0,0.0,6.0,6.0,0.0,0.0 };
  static double Tr[6] = { 0.0,0.0,0.0,10.0,0.0,0.0 };
  static double bem[6][6] =
   { { -1.072,-0.214,-1.072, ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 3.000,-1.072, 0.000, 0.000 },
     { -1.072, 3.000,-1.072, 1.286, 0.000, 0.000 },
     {  ATBOND,-1.072, 1.286,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static trna_loop rhit[NH];
  gene te,*tn;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,0,0,0,13,8,0,28,0,0,3,0,5,7,
     tmRNA,0.0,0,0,0 };
  tpos = th->pos;
  flag = 0;
  te.energy = sw->tmrnathresh;
  athresh = sw->tmathresh;
  cthresh = sw->tmcthresh;
  cathresh = sw->tmcathresh;
  s = tpos + tarm + 4;
  v = tpos + th->stem - 10;
  energy = K[*v] + G[v[1]] + A[v[2]];
  e = K[v[1]] + G[v[2]] + A[v[3]];
  if (e > energy) energy = e;
  if (energy < 18.0) energy = 0.0;
  tenergy = Tr[*s]+Cr[s[1]]+Cr[s[2]]+Ar[s[3]]+ energy + 1.59*the;
  na = -1;
  while (++na < nah)
   { aenergy = ahit[na].energy;
     if (aenergy < athresh) continue;
     apos = ahit[na].pos;
     if (apos < (tpos + MINTSTEM_DIST)) continue;
     if (apos > (tpos + MAXTSTEM_DIST + MAXPPINTRONDIST)) break;
     energy = -INACTIVE;
     sa = apos + t.astem1;
     for (sb=sa+9, se=sb+t.cstem; sb <= (sa+16); sb++,se++)
      for (sf = tpos-3; sf >= (tpos-7); sf--)
       { s1 = sb;
         s2 = sf;
         e = bem[*s1++][*--s2];
         while (s1 < se) e += bem[*s1++][*--s2];
         if (e > energy)
          { energy = e;
            t.var = (int)(tpos - sf);
            t.dloop = (int)(sb - sa); }}
     if (energy < cthresh) continue;
     energy += aenergy;
     if (energy < cathresh) continue;
     sb = sa + 3;
     sf = sa + 7;
     r = gtem[*sb++];
     while (sb < sf)
      { r = (r >> 4) + gtem[*sb++];
        if ((r & 3) == 2)
         { energy += 14.0;
           break; }}
     penergy = tenergy + Ga[apos[1]] + Ga[apos[2]] + energy;
     nrh = find_resume_seq(apos+MINTPDIST,TPWINDOW,rhit,NH,sw);
     nr = -1;
     while (++nr < nrh)
      { ps = rhit[nr].pos;
        t.energy = penergy + rhit[nr].energy;
        if (rhit[nr].stem < 24) t.energy -= 15.0;
        if (t.energy > te.energy)
         { flag = 1;
           t.tstem = th->stem;
           t.tloop = th->loop;
           t.asst = (long)(apos - tpos) + t.var + t.cstem;
           t.ps = tpos - t.var - t.cstem;
           t.tps = (int)(ps - t.ps);
           t.tpe = t.tps + rhit[nr].stem;
           te = t; }}}
  if (flag)
   { te.start = (long)(te.ps - seq) - 54;
     te.intron = te.cstem + te.var + 2*te.tstem + te.tloop +
                 te.astem2;
     as = aatail(te.ps + te.intron,&aext,sw);
     te.aatail = aext;
     base_copy3(te.ps-54,te.eseq,te.tpe+1+TMPTRAILER);
     te.nbase = te.astem1 + te.dloop + te.cstem;
     base_copy3(te.ps+te.asst,te.seq,te.nbase);
     base_copy3(te.ps,te.seq+te.nbase,te.intron + ASTEM2_EXT);
     te.intron += aext;
     te.nbase += te.intron;
     te.nintron = te.tpe - te.nbase + 1 + TMPTRAILER;
     te.intron += 54;
     te.tps += 54;
     te.tpe += 54;
     te.asst += 54;
     tn = find_slot(d,&te,&nts,sw);
     if (tn) *tn = te; }
  return(nts); }
  
 
 
  
  
int ti_genedetected(data_set *d, int nts, int *seq, gene *te, csw *sw)
{ int as,aext,as8,aext8,nbasefext,*s;
  int pseq[2*MAXETRNALEN+1];
  gene *tn;
  te->nbase = te->astem1 + te->spacer1 + te->spacer2 + 2*te->dstem +
              te->dloop +  2*te->cstem + te->cloop +
              te->var + 2*te->tstem + te->tloop + te->astem2;
  s = te->ps + te->nbase + te->nintron;
  as = aatail(s,&aext,sw);
  if (sw->extastem)
   if (te->astem1 == 7)
    if (bp[te->ps[-1]][*s])
     { as8 = aatail(s+1,&aext8,sw);
       if (as8 >= as)
        { te->ps--;
          te->nbase += 2;
          te->anticodon++;
          if (te->nintron > 0) te->intron++;
          te->astem1 = 8;
          te->astem2 = 8;
          as = as8;
          aext = aext8; }}
  nbasefext = te->nbase + ASTEM2_EXT;
  te->nbase += aext;
  te->start = (long)(te->ps - seq);
  tn = find_slot(d,te,&nts,sw);
  if (tn)
   { if (te->nintron == 0)
      base_copy3(te->ps,te->seq,nbasefext);
     else
      { base_copy3(te->ps,te->eseq,nbasefext + te->nintron);
        remove_intron(te->ps,pseq,nbasefext,
                      te->intron,te->nintron);
        base_copy3(pseq,te->seq,nbasefext); }
     te->aatail = aext;
     *tn = *te; }
  return(nts); }


int tmioptimise(data_set *d, int *seq, int lseq, int nts, csw *sw)
{ int i,j,k,intron,nt,nth,nd1,nd2,ndx,ndh,na,nah,nppah,nc,nch,tfold,tarm;
  int dstem,dloop,flag,mindist,maxdist,tmindist,tmaxdist,tmmindist,tmmaxdist;
  int tarmthresh,tmstrict,sp2min,sp2max,ige[7];
  int *se,*sc,*sb,*si,*tpos,*tend,*apos,*dpos,*tloopfold,*tmv,*cend;
  int *s1,*s2,*sd,*sf,*sl,*sg1,*sg2,*cposmin,*cposmax,*cpos;
  unsigned int r,q,c;
  double e,ec,he,the,thet,ethresh,energy,cenergy,denergy,ienergy;
  double tdarmthresh,genergy,energy2,energyf,energyf6;
  static unsigned int TT[6] =
   { 0x00, 0x00, 0x00, 0x11, 0x00, 0x00 };
  static unsigned int GG[6] =
   { 0x00, 0x00, 0x11, 0x00, 0x00, 0x00 };
  static unsigned int ct[6] = { 0,0,0,0,0,0 };
  static unsigned int cA[6] = { 0,0,0,2,0,0 };
  static unsigned int cC[6] = { 0,0,2,0,0,0 };
  static unsigned int cG[6] = { 0,2,0,1,0,0 };
  static unsigned int cT[6] = { 2,0,1,0,0,0 };
  static int yic[9] = { 1,0,0,0,0,0,0,0,0 };
  static int tic[9] = { 1,1,0,0,0,0,0,0,0 };
  static int a1ic[9] = { 1,1,1,0,0,0,0,0,0 };
  static int a2ic[9] = { 1,1,1,1,0,0,0,0,0 };
  static int a3ic[9] = { 1,1,1,1,1,0,0,0,0 };
  static int ric[9] = { 1,1,1,1,1,1,0,0,0 };
  static int goffb[13] = { 0,0,0,0,1,2,2,2,2,2,2,2,2 };
  static int goffe[13] = { 0,0,0,0,2,3,4,4,5,6,6,6,6 };
  static int cY[6] = { 0,1,0,1,0,0 };
  static int cR[6] = { 1,0,1,0,0,0 };
  static double ilw = 0.002;
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,7.0,0.0,0.0 };
  static double Y[6] = { 0.0,3.0,0.0,3.0,0.0,0.0 };
  static double R[6] = { 2.0,0.0,2.0,0.0,0.0,0.0 };
  static double YP[6] = { 0.0,3.0,0.0,3.0,0.0,0.0 };
  static double RP[6] = { 2.0,0.0,2.0,0.0,0.0,0.0 };
  static double RI[6] = { 0.1,0.0,0.05,0.0,0.0,0.0 };
  static double GI[6] = { 0.0,0.0,0.1,0.0,0.0,0.0 };
  static double YI[6] = { 0.0,0.1,0.0,0.1,0.0,0.0 };
  static double AI[6] = { 1.0,0.0,0.0,0.0,0.0,0.0 };
  static double GC[6] = { 0.0,1.5,6.0,0.0,0.0,0.0 };
  static double G3[6] = { 0.0,6.0,12.0,12.0,0.0,0.0 };
  static double dR[6] = { 6.0,0.0,6.0,0.0,0.0,0.0 };
  static double RH[6] = { 3.0,0.0,3.0,0.0,0.0,0.0 };
  static double AGT[6] = { 6.0,0.0,6.0,6.0,0.0,0.0 };
  static double dT[6] = { 0.0,0.0,0.0,6.0,0.0,0.0 };
  static double dbem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double dfem[6][6] =
   { { -4.000,-4.000,-4.000, ATBOND, 0.000, 0.000 },
     { -4.000,-4.000, 3.000,-4.000, 0.000, 0.000 },
     { -4.000, 3.000,-4.000, 1.286, 0.000, 0.000 },
     {  ATBOND,-4.000, 1.286,-4.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double cbem[6][6] =
   { { -1.072,-0.214,-1.072,2.0*ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 6.000,-1.072, 0.000, 0.000 },
     { -1.072, 6.000,-1.072, 3.400, 0.000, 0.000 },
     {  2.0*ATBOND,-1.072, 3.400,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };    
  static trna_loop thit[NTH],chit[NC],ahit[NA];
  static trna_dloop dhit[ND];
  gene te;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,2,1,3,9,5,7,0,0,0,15,0,5,7,
     tRNA,0.0,0,0,0 };  
  if (sw->mtrna)
   { nts = find_mt_trna(d,seq,lseq,nts,sw);
     if (!sw->tmrna) return(nts); }
  ethresh = sw->trnathresh;
  tmmindist = MINTPTSDIST + MINTPDIST;
  tmmaxdist = MAXTPTSDIST + MAXTPDIST;
  tmindist = (MINTRNALEN + sw->minintronlen - MAXTSTEM_DIST);
  tmaxdist = (MAXTRNALEN + sw->maxintronlen - MINTSTEM_DIST);
  if (sw->trna)
   { if (sw->tmrna)
      { mindist = (tmindist < tmmindist)?tmindist:tmmindist;
        maxdist = (tmaxdist > tmmaxdist)?tmaxdist:tmmaxdist; }
     else
      { mindist = tmindist;
        maxdist = tmaxdist; }}
  else
   { mindist = tmmindist;
     maxdist = tmmaxdist; }
  tarmthresh = sw->ttarmthresh;
  tdarmthresh = sw->tdarmthresh;
  tmstrict = sw->tmstrict;
  sp2min = sw->sp2min;
  sp2max = sw->sp2max;
  nth = find_tstems(seq,lseq,thit,NTH,sw);
  nt = -1;
  while (++nt < nth)
   { tpos = thit[nt].pos;
     t.tloop = thit[nt].loop;
     t.tstem = thit[nt].stem;
     tfold = tpos[-1];
     tloopfold = tpos + t.tstem + 1;
     tarm = 2*t.tstem + t.tloop;
     tend = tpos + tarm;
     tmv = tpos - VARMIN;
     flag = 0;
     te.energy = ethresh;
     the = thit[nt].energy;
     nah = find_astem5(tpos-maxdist,tpos-mindist,tend,7,ahit,NA,sw);
     if (sw->tmrna)
      { thet = the - G[tpos[t.tstem]] - G[tpos[t.tstem+1]];
        if (tmstrict)
         { if (thet >= tarmthresh)
            nts = tmopt(d,thit+nt,tarm,thet,ahit,nah,nts,seq,sw); }
        else
         nts = tmopt(d,thit+nt,tarm,the,ahit,nah,nts,seq,sw);
        nppah = find_astem5(tpos+MINPPASDIST,tpos+MAXPPASDIST,
                            tend,7,ahit+nah,NA-nah,sw);
        nts = tmopt_perm(d,thit+nt,tarm,the,ahit+nah,nppah,nts,seq,sw);
        if (thet < tarmthresh) continue;
        the = thet; }
     else if (sw->threshlevel < 1.0) /* find_tstems is generating extra tstems */
      { the -= (G[tpos[t.tstem]] + G[tpos[t.tstem+1]]);
        if (the < tarmthresh) continue; }   
     if (!sw->trna) continue;
     na = -1;
     while (++na < nah)
      { apos = ahit[na].pos;
        if (apos < (tpos - tmaxdist)) continue;
        if (apos > (tpos - tmindist)) break;
        he = the + ahit[na].energy;
        
  /* find dstems */
          
  ndh = 0;
  sc = apos + 8;
  energyf = dfem[sc[5]][tfold];
  sl = sc + sw->sp1max;
  while (sc < sl)
   { energy2 = dT[sc[-2]] + RH[*(sc-1)] + GC[*sc] + dfem[sc[-2]][sc[4]];
     energyf6 = dfem[sc[6]][tfold];
     for (dstem = 3; dstem <= 4; dstem++)
      { sd = sc + dstem;
        dloop = 3;
        se = sd + dloop;
        energy = energy2 + 6.0 + dR[*(se-1)] + energyf;
        if (dstem == 3)
         if (energyf < 0.0) energyf = energyf6;
        se += dstem;
        s1 = sc;
        s2 = se;
        sf = s1 + dstem;
        while (s1 < sf) energy += dbem[*s1++][*--s2];
        if (energy >= tdarmthresh)
         { if (ndh >= ND) goto DFL;
           dhit[ndh].pos = sc;
           dhit[ndh].end = se;
           dhit[ndh].loop = dloop;
           dhit[ndh].stem = dstem;
           dhit[ndh].energy = energy;
           ndh++; }
        sg1 = sd + 1;
        sg2 = sd + 6;
        q = GG[*sg1++];
        ige[1] = q & 3;
        j = 2;
        while (sg1 <= sg2)
         { q = (q >> 4) + GG[*sg1++];
           ige[j++] = q & 3; }
        for (dloop = 4; dloop <= 11; dloop++)
         { j = goffb[dloop];
           k = goffe[dloop];
           c = ige[j++];
           while (j <= k) c = c | ige[j++];
           genergy = G3[c];
           se = sd + dloop;
           energy = energy2 + genergy + dR[*(se-1)] + energyf;
           se += dstem;
           s1 = sc;
           s2 = se;
           sf = s1 + dstem;
           while (s1 < sf) energy += dbem[*s1++][*--s2];
           if (energy >= tdarmthresh)
            { if (ndh >= ND) goto DFL;
              dhit[ndh].pos = sc;
              dhit[ndh].end = se;
              dhit[ndh].loop = dloop;
              dhit[ndh].stem = dstem;
              dhit[ndh].energy = energy;
              ndh++; }}}
      s1 = sc;
      s2 = sc + 16;
      sd = sc + 6;
      j = bp[*s1][*--s2];
      while (++s1 < sd) j += bp[*s1][*--s2];
      if (j >= 6)
       { energy = dT[sc[-1]] + RH[*sc] + GC[*(sc+1)] + energyf6;
         energy += G[*++sd];
         energy += G[*++sd];
         energy += AGT[*++sd] + dfem[sc[-1]][sc[4]];
         sd += 7;
         s1 = sc;
         s2 = sd;
         sf = s1 + 6;
         while (s1 < sf) energy += dbem[*s1++][*--s2];
         if (energy >= tdarmthresh)
          { if (ndh >= ND) goto DFL;
            dhit[ndh].pos = sc;
            dhit[ndh].end = sd;
            dhit[ndh].loop = 4;
            dhit[ndh].stem = 6;
            dhit[ndh].energy = energy;
            ndh++; }}
      s1 = sc;
      s2 = sc + 18;
      sd = sc + 7;
      j = bp[*s1][*--s2];
      while (++s1 < sd) j += bp[*s1][*--s2];
      if (j >= 7)
       { energy = energy2 + dfem[sc[7]][tfold];
         energy += G[*++sd];
         energy += G[*++sd];
         energy += AGT[*++sd];
         sd += 8;
         s1 = sc;
         s2 = sd;
         sf = s1 + 7;
         while (s1 < sf) energy += dbem[*s1++][*--s2];
         if (energy >= tdarmthresh)
          { if (ndh >= ND) goto DFL;
            dhit[ndh].pos = sc;
            dhit[ndh].end = sd;
            dhit[ndh].loop = 4;
            dhit[ndh].stem = 7;
            dhit[ndh].energy = energy;
            ndh++; }}
     energyf = energyf6;
     sc++; }
  goto DFN;
  DFL:
  fprintf(stderr,"Too many D-stem hits\n");
  DFN:

  /* End of find dstems routine */                       
                          
        nd1 = ndh;
        while (--nd1 >= 0)
         { dstem = dhit[nd1].stem;
           dpos = dhit[nd1].pos;          
           if ((int)(dpos - apos) < 9) 
             dhit[nd1].energy -= 3.0;          
           if (*tloopfold == Guanine)
            { sb = dpos + dstem + 2;
              sc = sb;
              se = sb + dhit[nd1].loop - 3;
              r = TT[*sb++];
              while (sb < se)
               { r = (r >> 4) + TT[*sb++];
                 if (r & 2)
                  { dhit[nd1].energy += 10.0;
                    break; }}
              r = GG[*sc++];
              while (sc < se)
               { r = (r >> 4) + GG[*sc++];
                 if (r & 2)
                 { dhit[nd1].energy -= 12.0;
                   break; }}}}
        nd1 = ndh;
        while (--nd1 >= 0)
         { if (!dhit[nd1].end) continue;
           cpos = dhit[nd1].end;
           denergy = dhit[nd1].energy;
           ndx = nd1;
           nd2 = nd1;
           while (--nd2 >= 0)
            { if (dhit[nd2].end != cpos) continue;
              e = dhit[nd2].energy;
              if (e > denergy)
               { denergy = e;
                 dhit[ndx].end = NULL;
                 ndx = nd2; }}}
        cposmin = 0;
        cposmax = 0;
        nd1 = ndh;
        while (--nd1 >= 0)
         { if (!dhit[nd1].end) continue;
           cposmin = dhit[nd1].end;
           cposmax = cposmin;
	       break; }
        nd2 = nd1;
        while (--nd2 >= 0)
         { if (!(cpos = dhit[nd2].end)) continue;
           if (cpos < cposmin) cposmin = cpos;
           if (cpos > cposmax) cposmax = cpos; }
	for (cpos = cposmin + sp2min; cpos <= (cposmax + sp2max); cpos++)
         { denergy = -INACTIVE;
           ndx = -1;
           nd1 = ndh;
           while (--nd1 >= 0)
            { if (!dhit[nd1].end) continue;
              if ((dhit[nd1].end + sp2max) < cpos) continue; 
              if ((dhit[nd1].end + sp2min) > cpos) continue; 
              e = dhit[nd1].energy;
              if (e > denergy)
               { denergy = e;
                 ndx = nd1; }}
           if (ndx < 0) continue;
           denergy += he;
           if (denergy < (te.energy - 49.0)) continue;
           
  /* find cstems */

  nch = 0;
  si = cpos;
  sc = cpos + 5;
  se = cpos + 4;
  ct[0] = cA[*se];
  ct[1] = cC[*se];
  ct[2] = cG[*se];
  ct[3] = cT[*se];
  while (--se >= cpos)
   { ct[0] = (ct[0] << 4) + cA[*se];
     ct[1] = (ct[1] << 4) + cC[*se];
     ct[2] = (ct[2] << 4) + cG[*se];
     ct[3] = (ct[3] << 4) + cT[*se]; }
  si += 11;
  se = tmv - VARDIFF - 5;
  if (si < se) si = se;
  r = ct[*si++];
  r = (r >> 4) + ct[*si++];
  r = (r >> 4) + ct[*si++];
  r = (r >> 4) + ct[*si++];
  while (si < tmv)
   { r = (r >> 4) + ct[*si++];
     if ((r & 0xf) >= 5)
      { if (nch >= NC)
         { fprintf(stderr,"Too many cstem hits\n");
           goto FN; }
        chit[nch].pos = si;
        chit[nch].stem = 5;
        chit[nch].loop = (int)(si - sc - 5);
        if (chit[nch].loop == 9)
         if (bp[*sc][si[-6]])
          if (cY[sc[2]])
           if (cR[sc[6]])
            if (cY[sc[1]])
             { chit[nch].stem = 6;
               chit[nch].loop = 7; }
        s1 = cpos;
        s2 = si;
        se = s1 + chit[nch].stem;
        chit[nch].energy = cbem[*s1++][*--s2];
        while (s1  < se)
         chit[nch].energy += cbem[*s1++][*--s2];
        nch++; }}
  FN:

  /* end of find cstems routine */

           nc = -1;
           while (++nc < nch)
            { energy = denergy + chit[nc].energy;
              if (energy < (te.energy - 19.0)) continue;
              cend = chit[nc].pos;
              t.var = (int)(tpos - cend);
              t.cloop = chit[nc].loop;
              t.cstem = chit[nc].stem;
              intron = 0;
              if (t.cloop < 9)
               { if (sw->minintronlen > 0) continue;
		         if (sw->cloop7)
		          if (t.cloop != 7) continue;
                 t.nintron = 0;
                 if (t.var > 17) energy += vloop_stability(cend,t.var,&t.varbp);
                 sb = cpos + t.cstem;
                 energy += T[*(sb + 1)] + Y[*(sb)] + R[*(sb + 5)] -
                           0.05*t.var - ((t.cloop == 7)?0.0:6.0); }
              else
               { t.nintron = t.cloop - 7;
                 if (t.nintron > sw->maxintronlen) continue;
                 if (t.nintron < sw->minintronlen) continue;
                 if (t.var > 17) energy += vloop_stability(cend,t.var,&t.varbp);
                 if (energy < (te.energy - 9.0)) continue;
                 t.cloop = 7;
                 sb = cpos + t.cstem;
                 se = sb + t.nintron;
		 if (sw->ifixedpos)
	          { intron = 6;
                    cenergy = YP[*sb] + T[sb[1]] + RP[sb[5]]; }
                 else
                  { cenergy = YP[*se] + T[*(se+1)] + RP[*(se+5)];
                    ienergy = cenergy + RI[*sb] + GI[*(se-1)] +
                              AI[se[-2]]*YI[se[-1]];
                    for (j = 1; j <= 7; j++)
                     { si = se + j - 1;
                       ec = YP[*(sb + yic[j]*t.nintron)] +
                           T[*(sb + tic[j]*t.nintron + 1)] +
                           RP[*(sb + ric[j]*t.nintron + 5)];
                       e = ec + RI[*(sb + j)] + GI[*si] + AI[si[-1]]*YI[*si];
                       if (j == 6) e += 0.01;
                       if (e > ienergy)
                        { ienergy = e;
                          cenergy = ec;
                          intron = j; }}}
                  energy +=  cenergy - 10.0 - ilw*(t.nintron  + 1.1*t.var);
                  if (t.nintron >= 130)
                   { si = se + intron;
                     j = si[-1];
                     if (j != Guanine)
                      { if (si[-2] != Adenine) energy -= 4.0;
                        if (j != Cytosine)
                         if (j != Thymine)
                          energy -= 8.0; }}}
              dstem = dhit[ndx].stem;
              dpos = dhit[ndx].pos;
              if (dstem >= 6)
               { if (sb[2 + a1ic[intron]*t.nintron] != Thymine) continue;
                 if (sb[3 + a2ic[intron]*t.nintron] != Cytosine) continue;
                 if (sb[4 + a3ic[intron]*t.nintron] != Adenine) continue;
                 energy += 3.0; }
              else
               if (!(dpos[-1] & 5))
                { i = 0;
                  si = cend;
                  se = cend + 4;
                  while (si < se)
                   { if (!(*si++ & 5))
                      { if (++i >= 2)
                         { energy += 3.0;
                           break; }}
                     else
                      i = 0; }}
              if (t.cstem >= 6)
               { if (sb[2 + a1ic[intron]*t.nintron] == Cytosine)
                  if (sb[3 + a2ic[intron]*t.nintron] == Thymine)
                   if (sb[4 + a3ic[intron]*t.nintron] == Adenine) 
                    energy += 4.0; }
          if (energy < ethresh) continue;
          t.energy = energy;
          t.dstem = dstem;
          t.astem1 = (t.dstem < 6)?7:((t.tstem < 5)?9:8);
          t.astem2 = t.astem1;
          t.ps = apos + 7 - t.astem1;
          t.nbase = (int)(tend - t.ps) + t.astem2;
          t.dloop = dhit[ndx].loop;
          t.spacer1 = (int)(dpos - apos) - 7;
          t.spacer2 = (int)(cpos - dhit[ndx].end);
          j = (int)(cpos - t.ps) + t.cstem;
          t.anticodon = j + 2;
          if (t.nintron > 0)
           { t.intron = j + intron;
             if ((t.nbase + t.nintron) > MAXTRNALEN)
              { nts = ti_genedetected(d,nts,seq,&t,sw);
                continue; }}
          if (energy < te.energy) continue;
          flag = 1;
          te = t; } }}
     if (flag) nts = ti_genedetected(d,nts,seq,&te,sw); }
  return(nts); }


void disp_ftable_entry(FILE *f, int n[], int i, int m, csw *sw)
 { if (m > 0)
    switch(sw->geneticcode)
      { case METAZOAN_MT:
          fprintf(f," %-18s %-4d",aa(n,sw),m);
          break;
        case STANDARD:
        case VERTEBRATE_MT:
        default:
          fprintf(f," %-4s %-5d",aa(n,sw),m);
          break; }
   else
    switch(sw->geneticcode)
     { case METAZOAN_MT:
         fprintf(f," %-18s     ",aa(n,sw));
         break;
      case STANDARD:
      case VERTEBRATE_MT:
      default:
        fprintf(f," %-4s      ",aa(n,sw));
        break; }}


void disp_freq_table(int nt, csw *sw)
{ int i,j,k,m,ambig,*s,c1,c2,c3,c[3],a[3],table[4][4][4];
  static int cgflip[4] = { 0,2,1,3 };
  static int codonorder[4] = { 3,1,0,2 };
  FILE *f = sw->f;
  ambig = 0;
  for (i = 0; i < 4; i++)
   for (j = 0; j < 4; j++)
    for (k = 0; k < 4; k++)
     table[i][j][k] = 0;
  for (i = 0; i < nt; i++)
   if (ts[i].energy >= 0.0)
    if (ts[i].genetype == tRNA)
     if (ts[i].cloop == 7)
     { s = ts[i].seq + ts[i].anticodon;
       c1 = *s;
       c2 = s[1];
       c3 = s[2];
       if ((c1 >= Adenine) && (c1 <= Thymine))
        if ((c2 >= Adenine) && (c2 <= Thymine))
         if ((c3 >= Adenine) && (c3 <= Thymine))
          table[*s][s[1]][s[2]]++;
         else ambig++;
        else ambig++;
       else ambig++; }
     else ambig++;
  fprintf(f,"tRNA anticodon frequency\n");
  for (i = 0; i < 4; i++)
   { c[0] = codonorder[i];
     a[2] = 3 - c[0];
     for (j = 0; j < 4; j++)
      { c[2] = codonorder[j];
        a[0] = 3 - c[2];
        for (k = 0; k < 4; k++)
         { c[1] = codonorder[k];
           a[1] = 3 - c[1];
           fprintf(f,"%c%c%c",cpbase(a[0]),cpbase(a[1]),cpbase(a[2]));
           m = table[a[0]][a[1]][a[2]];
           disp_ftable_entry(f,a,k,m,sw); }
        fputc('\n',f); }
     if (i < 3) fputc('\n',f); }
  if (ambig > 0) fprintf(f,"Ambiguous: %d\n",ambig);
  fprintf(f,"\ntRNA codon frequency\n");
  for (i = 0; i < 4; i++)
   { c[0] = codonorder[i];
     a[2] = 3 - c[0];
     for (j = 0; j < 4; j++)
      { c[2] = codonorder[j];
        a[0] = 3 - c[2];
        for (k = 0; k < 4; k++)
         { c[1] = codonorder[k];
           a[1] = 3 - c[1];
           fprintf(f,"%c%c%c",cpbase(c[0]),cpbase(c[1]),cpbase(c[2]));
           m = table[a[0]][a[1]][a[2]];
           disp_ftable_entry(f,a,k,m,sw); }
        fputc('\n',f); }
     if (i < 3) fputc('\n',f); }
  if (ambig > 0) fprintf(f,"Ambiguous: %d\n",ambig);
  fputc('\n',f); }

void disp_energy_stats(data_set *d, int nt, csw *sw)
{ int i,n[NS],genetype,introns,nintron,trna,mtrna,ntv,nd,nps;
  double gc,gcmin[NS],gcmax[NS];
  FILE *f = sw->f;
  mtrna = sw->mtrna;
  trna = sw->trna | mtrna;
  nps = 0;
  if (mtrna)
   { ntv = 0;
     nd = 0; }
  if ((sw->trna) && (sw->maxintronlen > 0))
   { introns = 1;
     nintron = 0; }
  else
   introns = 0;
  for (i = 0; i < NS; i++)
   { n[i] = 0;
     gcmin[i] = 1.0;
     gcmax[i] = 0.0; }
  for (i = 0; i < nt; i++)
   if (ts[i].energy >= 0.0)
    { n[NS-1]++;
      genetype = ts[i].genetype;
      n[genetype]++;
      if (pseudogene(ts + i,sw)) nps++;
      if (genetype == tRNA)
       { if (mtrna)
          { if (ts[i].tstem == 0) ntv++;
            if (ts[i].dstem == 0) nd++; }
         if (introns) if (ts[i].nintron > 0) nintron++;
         gc = gc_content(ts+i);
         if (gc < gcmin[genetype]) gcmin[genetype] = gc;
         if (gc > gcmax[genetype]) gcmax[genetype] = gc; }}
  fputc('\n',f);
  fputc('\n',f);
  if (sw->repeatsn)
   if ((n[tRNA] + n[tmRNA]) > 0)
     fprintf(f,"%s\n\n",d->seqname);
  if (trna)
   { sw->ngene[tRNA] += n[tRNA];
     if (n[tRNA] > 3) disp_freq_table(nt,sw);
     if ((n[tRNA] > 1) || ((sw->tmrna) && (n[tmRNA] > 0)))
      { if (introns)
         { if (sw->minintronlen == 0)
            fprintf(f,"Number of tRNA genes with no introns = %d\n",
                    n[0]-nintron);
           fprintf(f,"Number of tRNA genes with C-loop introns = %d\n",
                   nintron); }
        else
         fprintf(f,"Number of %s genes = %d\n",sw->genetypename[tRNA],n[tRNA]);
        if (mtrna)
         { if (sw->tvloop)
            fprintf(f,"Number of TV replacement loop tRNA genes = %d\n",
             ntv);
           fprintf(f,"Number of D replacement loop tRNA genes = %d\n",
     nd); }
  if (n[tRNA] > 1)
   fprintf(f,"tRNA GC range = %2.1f%% to %2.1f%%\n",
           gcmin[0]*100.0,gcmax[0]*100.0); }}
  if (sw->tmrna)
   { sw->ngene[tmRNA] += n[tmRNA];
     if ((n[tmRNA] > 1) || (trna && (n[tRNA] > 0)))
      fprintf(f,"Number of %s genes = %d\n",sw->genetypename[tmRNA],n[tmRNA]); }
  sw->nps += nps;
  if (sw->reportpseudogenes)
   if (nps > 0) 
    if (n[NS-1] > 1)
     fprintf(f,"Number of possible pseudogenes = %d\n",nps);
  fputc('\n',f);
  fputc('\n',f); }
  

void batch_energy_stats(data_set *d, int nt, csw *sw)
{ int i,n[NS],genetype,introns,nintron,trna,mtrna,ntv,nd,nps;
  double gc,gcmin[NS],gcmax[NS];
  FILE *f = sw->f;
  mtrna = sw->mtrna;
  trna = sw->trna | mtrna;
  nps = 0;
  if (mtrna)
   { ntv = 0;
     nd = 0; }
  if ((sw->trna) && (sw->maxintronlen > 0))
   { introns = 1;
     nintron = 0; }
  else
   introns = 0;
  for (i = 0; i < NS; i++)
   { n[i] = 0;
     gcmin[i] = 1.0;
     gcmax[i] = 0.0; }
  for (i = 0; i < nt; i++)
   if (ts[i].energy >= 0.0)
    { n[NS-1]++;
      genetype = ts[i].genetype;
      n[genetype]++;
      if (ts[i].energy < 100.0) nps++;
      if (genetype == tRNA)
       { if (mtrna)
          { if (ts[i].tstem == 0) ntv++;
            if (ts[i].dstem == 0) nd++; }
         if (introns) if (ts[i].nintron > 0) nintron++;
         gc = gc_content(ts+i);
         if (gc < gcmin[genetype]) gcmin[genetype] = gc;
         if (gc > gcmax[genetype]) gcmax[genetype] = gc; }}
  if (trna) sw->ngene[tRNA] += n[tRNA];
  if (sw->tmrna) sw->ngene[tmRNA] += n[tmRNA];
  sw->nps += nps; }


int gene_sort(data_set *d, int nt, int sort[], csw *sw)
{ int i,n,j,k;
  long starti,startj,stopi,stopj,psmax;
  psmax = d->psmax;
  n = 0;
  for (i = 0; i < nt; i++)
   if (ts[i].energy >= 0.0)
    { if (sw->ireportminintronlen == 1)
       if (ts[i].genetype == tRNA)
        if (ts[i].nintron < sw->minintronlenreport)
         continue;
      sort[n++] = i; }
  i = -1;
  while (++i < (n-1))
   { j = i;
     while (++j < n)
      { starti = ts[sort[i]].start;
        startj = ts[sort[j]].start;
        stopi = ts[sort[i]].stop;
        stopj = ts[sort[j]].stop;
        if (stopi < starti)
         if ((psmax - starti) < stopi) starti -= psmax;
         else stopi += psmax;
        if (stopj < startj)
         if ((psmax - startj) < stopj) startj -= psmax;
         else stopj += psmax;
        if (starti > startj)
         { k = sort[i];
           sort[i] = sort[j];
           sort[j] = k; }
        else
         if (starti == startj)
          if (stopi < stopj)
           { k = sort[i];
             sort[i] = sort[j];
             sort[j] = k; }}}
   return(n); }


int iamatch(data_set *d, gene *t, csw *sw)
{ char key[5],*k,s[100];
  if (k = softstrpos(d->seqname,"TRNA-")) k += 5;
  else
   if (k = wildstrpos(d->seqname,"|***|")) k++;
    else return(-1);
  copy3cr(k,key,3);
  name(t,s,1,sw);
  if (softstrpos(s,key)) return(1);
  return(0); }



int gene_mismatch(data_set *d, annotated_gene *agene, gene *t, csw *sw)
{
int w,alen,dlen;
char *s;
w = 0;
dlen = seqlen(t);
alen = aseqlen(d,agene);
switch(t->genetype)
 { case tRNA:
     s = aa(t->seq + t->anticodon,sw);
     if (!softstrpos(s,agene->species+5))
      { if (t->cloop == 8)
         { s = aa(t->seq + t->anticodon + 1,sw);
           if (!softstrpos(s,agene->species+5)) w += 1; }
        else if (t->cloop == 6)
         { s = aa(t->seq + t->anticodon - 1,sw);
           if (!softstrpos(s,agene->species+5)) w += 1; }
        else w += 1; }
     if (agene->comp != t->comp) w += 2;
     if (alen <= (dlen - sw->trnalenmisthresh)) w += 4;
     else if (alen >= (dlen + sw->trnalenmisthresh)) w += 4;
     break;
  case tmRNA:
     if (agene->comp != t->comp) w += 2;
     if (alen <= (dlen - sw->tmrnalenmisthresh)) w += 4;
     else if (alen >= (dlen + sw->tmrnalenmisthresh)) w += 4;
     break; }
return(w);
}


int gene_mismatch_report(data_set *d, annotated_gene *agene, gene *t, char *report, csw *sw)
{
int w;
char *s;
w = gene_mismatch(d,agene,t,sw);
s = report;
if (w & 1) s = copy("amino acceptor",s);
if (w & 2)
 { if (w & 1) 
    if (w & 4) s = copy(", ",s);
    else s = copy(" and ",s);
   s = copy("sense",s); }
if (w & 4) 
 { if ((w & 3) > 0) s = copy(" and ",s);
   s = copy("sequence length",s); }
if (w > 0) s = copy(" mismatch",s);
*s = '\0';
return(w);
}



int nearest_annotated_gene(data_set *d, gene *t,
                           int list[], int score[], int nmax,
                           csw *sw)
{ int n,i,j,k,q,w,nagene;
  long a,b,c,e,thresh,psmax;
  char *s; 
  annotated_gene *ta;
  psmax = d->psmax;
  nagene = d->nagene[NS-1];
  ta = d->gene;
  n = 0;
  a = t->start;
  b = t->stop;
  thresh = b-a;
  if (b < a)
   { b += psmax;
     thresh += psmax;
     for (i = 0; i < nagene; i++)
      { c = ta[i].start;
        e = ta[i].stop;
        if (e < c)
         { e += psmax;
           if (a > e) goto NXTW;
           if (b < c) goto NXTW;
           if (n >= nmax) break;
           list[n] = i;
           score[n] = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
           n++;
           NXTW:
           c -= psmax;
           e -= psmax; }
        if (a > e) continue;
        if (b < c) continue;
        if (n >= nmax) break;
        list[n] = i;
        score[n] = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
        n++; }
     a -= psmax;
     b -= psmax; }
  for (i = 0; i < nagene; i++)
   { c = ta[i].start;
     e = ta[i].stop;
     if (e < c)
      { e += psmax;
        if (a > e) goto NXTN;
        if (b < c) goto NXTN;
        if (n >= nmax) break;
        list[n] = i;
        score[n] = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
        n++;
        NXTN:
        c -= psmax;
        e -= psmax; }
     if (a > e) continue;
     if (b < c) continue;
     if (n >= nmax) break;
     list[n] = i;
     score[n] = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
     n++; }
  for (i = 0; i < n; i++)
   { k = list[i];
     if (ta[k].genetype == t->genetype) 
      { score[i] += 5000;
        w = gene_mismatch(d,ta + k,t,sw);
        if (w & 1) score[i] -= 2;
        if (w & 2) score[i] -= 1; }}
  if (n > 1)
   { for (i = 0; i < (n-1); i++)
      for (j = i+1; j < n; j++)
       if (score[j] > score[i])
       { k = list[i];
         list[i] = list[j];
         list[j] = k;
         k = score[i];
         score[i] = score[j];
         score[j] = k; }}
  return(n); }

    


int proximity_compare(data_set *d, int is,
                      long prox, long dlen, long alen,
                      annotated_gene *a,
                      csw *sw)
{
int w,score;
long diff;
char nm[200];
gene *t;
t = ts + is;
w = gene_mismatch(d,a,t,sw);
if (prox >= alen)
 { diff = dlen - alen;
   if (prox >= (2L*diff)) score = (int)(prox - diff);
   else score = (int)(prox/2L); }
else
 if (prox >= dlen)
  { diff = alen - dlen;
    if (prox >= (2L*diff)) score = (int)(prox - diff);
    else score = (int)(prox/2L); }
else { score = (int)prox; }
if (w & 1) score -= 10;
if (w & 2) score -= 2;
if (score < 0) score = 0;
if (t->annotation >= 0)
 if (t->annosc >= score) return(-1);
return(score);
}




int nearest_detected_gene(data_set *d, int sort[], int nd,
                          int *scorep,
                          annotated_gene *ag, csw *sw)
{ int n,i,is;
  long a,b,c,e,score,alen,scoremax,psmax;
  long prox,proximity;
  psmax = d->psmax;
  n = -1;
  scoremax = -1;
  a = ag->start;
  b = ag->stop;
  alen = b - a;
  if (b < a) alen += psmax;
  proximity = 1 + alen/2;
  if (proximity > 30) proximity = 30;
  if (b < a)
   { b += psmax;
     for (i = 0; i < nd; i++)
      { is = sort[i];
        if (ag->genetype != ts[is].genetype) continue;
        c = ts[is].start;
        e = ts[is].stop;
        if (e < c)
         { e += psmax;
           if (a > e) goto NXTW;
           if (b < c) goto NXTW;
           prox = (a >= c)?((b >= e)?e-a:alen):((b >= e)?e-c:b-c);
           if (prox >= proximity)
            if ((score = proximity_compare(d,is,prox,e-c,alen,ag,sw)) > scoremax)
             { n = i;
               scoremax = score; }
           NXTW:
           c -= psmax;
           e -= psmax; }
        if (a > e) continue;
        if (b < c) continue;
        prox = (a >= c)?((b >= e)?e-a:alen):((b >= e)?e-c:b-c);
        if (prox >= proximity)
         if ((score = proximity_compare(d,is,prox,e-c,alen,ag,sw)) > scoremax)
          { n = i;
            scoremax = score; }}
     a -= psmax;
     b -= psmax; }
  for (i = 0; i < nd; i++)
   { is = sort[i];
     if (ag->genetype != ts[is].genetype) continue;
     c = ts[is].start;
     e = ts[is].stop;
     if (e < c)
      { e += psmax;
        if (a > e) goto NXTN;
        if (b < c) goto NXTN;
        prox = (a >= c)?((b >= e)?e-a:alen):((b >= e)?e-c:b-c);
        if (prox >= proximity)
         if ((score = proximity_compare(d,is,prox,e-c,alen,ag,sw)) > scoremax)
          { n = is;
            scoremax = score; }
        NXTN:
        c -= psmax;
        e -= psmax; }
     if (a > e) continue;
     if (b < c) continue;
     prox = (a >= c)?((b >= e)?e-a:alen):((b >= e)?e-c:b-c);
     if (prox >= proximity)
      if ((score = proximity_compare(d,is,prox,e-c,alen,ag,sw)) > scoremax)
       { n = is;
         scoremax = score; }}
  *scorep = scoremax;
  return(n); }



void disp_match(data_set *d, int *sort, int nd, csw *sw)
{
int i,ld,fn[NS],fp[NS],fpd,fptv,w,score,detect,n[NS];
int prevannoted,nl,k,csort[NGFT],*msort;
long start;
char tag[52],nm[100],anm[100],ps[100],mreport[100],*s;
FILE *f = sw->f;
gene *t;
annotated_gene *agene,*a;
static char gp[2][7] = { "genes","gene" };
static char comp[3] = " c";
static char aps[2][5] = { "  ","PS" };
nl = nd;
if (sw->trna | sw->mtrna) nl += d->nagene[tRNA];
if (sw->tmrna) nl += d->nagene[tmRNA];
if (nl < NGFT) msort = csort;
else
 { msort = (int *)malloc(nl*sizeof(int));
   if (msort == NULL)
    { fprintf(stderr,"Not enough memory to match genes\n");
      return; }}
fprintf(f,"\n%s\n",d->seqname);
fprintf(f,"%ld nucleotides in sequence\n",d->psmax);
fprintf(f,"Mean G+C content = %2.1f%%\n",100.0*d->gc);
fprintf(f,"\nGenBank to Aragorn comparison\n\n");
sw->dispmatch = 1;
for (i = 0; i < NS; i++) 
  { n[i] = 0;
    fn[i] = 0;
    fp[i] = 0; }
for (i = 0; i < nd; i++)
 { w = sort[i];
   if (ts[w].energy >= 0.0)
    { n[NS-1]++;
     n[ts[w].genetype]++; }
   ts[w].annotation = -1;
   ts[w].annosc = -1; }
if (sw->trna | sw->mtrna | sw->tmrna)
 { fpd = 0;
   fptv = 0;
   if (sw->trna | sw->mtrna)
    { fprintf(f,"%d annotated tRNA %s\n",d->nagene[tRNA],gp[(d->nagene[tRNA]==1)?1:0]);
      fprintf(f,"%d detected tRNA %s\n",n[tRNA],gp[(n[tRNA]==1)?1:0]); }
   if (sw->tmrna)
    { fprintf(f,"%d annotated tmRNA %s\n",d->nagene[tmRNA],gp[(d->nagene[tmRNA]==1)?1:0]);
      fprintf(f,"%d detected tmRNA %s\n",n[tmRNA],gp[(n[tmRNA]==1)?1:0]); }
   fprintf(f,"\n  GenBank                                      Aragorn\n");
   nl = 0;
   for (i = 0; i < d->nagene[NS-1]; i++)
    { agene = d->gene + i;
      agene->detected = -1;
      if (agene->genetype != tRNA) 
       { if (agene->genetype != tmRNA) continue;
         else if (!sw->tmrna) continue; }
      else if (!sw->trna) if (!sw->mtrna) continue;
      a = agene;
      k = i;
      while ((a->detected = nearest_detected_gene(d,sort,nd,&score,a,sw)) >= 0)
       { t = ts + a->detected;
         prevannoted = t->annotation;
         t->annotation = k;
         t->annosc = score;
         if (prevannoted < 0) break;
         if (prevannoted == k) break;
         if (prevannoted == i) break; 
         a = d->gene + prevannoted;
         k = prevannoted; }
      k = nl;
      while (--k >= 0)
       { if (agene->start >= d->gene[msort[k]].start) break;
         msort[k+1] = msort[k]; }
      msort[++k] = i;
      nl++; }
   for (i = 0; i < nd; i++)
    { t = ts + sort[i];
      if (t->annotation >= 0) continue;
      if (t->genetype != tRNA) 
       { if (t->genetype != tmRNA) continue;
         else if (!sw->tmrna) continue; }
      else if (!sw->trna) if (!sw->mtrna) continue;
      k = nl;
      while (--k >= 0)
       { if (msort[k] >= 0) start = d->gene[msort[k]].start;
         else start = ts[-1-msort[k]].start;
         if (t->start >= start) break;
         msort[k+1] = msort[k]; }
      msort[++k] = -(sort[i] + 1);
      nl++; }
   for (i = 0; i < nl; i++)
    { if (msort[i] >= 0)
       { agene = d->gene + msort[i];
         detect = agene->detected;
         if (detect >= 0)
          { t = ts + detect;
            w = gene_mismatch_report(d,agene,t,mreport,sw);
            if (w > 0) fputc('*',f);
            else fputc(' ',f); }
         else fputc('*',f);
         sprintf(anm," %-11s%c(%ld,%ld) %s",
                 agene->species,comp[agene->comp],
                 sq(agene->start),sq(agene->stop),aps[agene->pseudogene]);
         fprintf(f,"%-45s ",anm);
         if (detect >= 0)
          { fprintf(f,"%s ",name(t,nm,1,sw));
            if (t->comp == 0) fputc(' ',f);
            fprintf(f,"%s",position(ps,t,sw));
            if (sw->energydisp) fprintf(f," %7.3lf",t->energy);
            if (t->genetype == tmRNA)
             { peptide_tag(tag,50,t,sw);
               fprintf(f," %s",tag); }
            if (sw->reportpseudogenes)
             if (pseudogene(t,sw))
              fprintf(f," PS");
            if (w > 0) fprintf(f," %s",mreport);
            fputc('\n',f); }
         else
          { fprintf(f,"Not detected\n");
            fn[agene->genetype]++; }}
      else
       { t = ts - (msort[i] + 1);
         fprintf(f,"* Not annotated                                %s ",name(t,nm,1,sw));
         if (t->comp == 0) fputc(' ',f);
         fprintf(f,"%s",position(ps,t,sw));
         if (sw->energydisp) fprintf(f," %7.3lf",t->energy);
         if (t->genetype == tmRNA)
          { peptide_tag(tag,50,t,sw);
            fprintf(f," %s",tag); }
         if (sw->reportpseudogenes)
          if (pseudogene(t,sw))
           fprintf(f," PS");
         fputc('\n',f);
         fp[t->genetype]++;
         if (t->genetype == tRNA)
          { if (t->dstem == 0) fpd++;
            if (t->tstem == 0) fptv++; }}}
   fprintf(f,"\n");
   if (sw->trna | sw->mtrna)
    { fprintf(f,"Number of annotated tRNA genes not detected = %d\n",fn[tRNA]);
      fprintf(f,"Number of unannotated tRNA genes detected = %d\n",fp[tRNA]); }
   if (sw->mtrna)
    { fprintf(f,"Number of unannotated D-replacement tRNA genes detected = %d\n",fpd);
      fprintf(f,"Number of unannotated TV-replacement tRNA genes detected = %d\n",fptv); }
   if (sw->tmrna)
    { fprintf(f,"Number of annotated tmRNA genes not detected = %d\n",fn[tmRNA]);
      fprintf(f,"Number of unannotated tmRNA genes detected = %d\n",fp[tmRNA]); }
   fprintf(f,"\n\n");
   for (i = tRNA; i <= tmRNA; i++)
    { sw->nagene[i] += d->nagene[i];
      sw->nafn[i] += fn[i]; 
      sw->nafp[i] += fp[i]; }
   if (sw->mtrna)
    { sw->natfpd += fpd;
      sw->natfptv += fptv; }}
sw->nabase += d->psmax;
sw->dispmatch = 0;
if (nl >= NGFT) free((void *)msort);
}


void annotation_overlap_check(data_set *d, gene *t, char *margin, csw *sw)
{
int a,m,n,w,list[20],score[20];
char mreport[100];
static char comp[3] = " c";
n = nearest_annotated_gene(d,t,list,score,20,sw);
if (n < 1) m = -1;
else
 { m = 0;
   a = list[m];
   if (d->gene[a].genetype != t->genetype) m = -1;
   else 
    { w = gene_mismatch_report(d,d->gene+a,t,mreport,sw);
      if (w & 1)
       { if ((score[m] - 5000) < (3*seqlen(t)/4)) m = -1; }
      else
       { if ((score[m] - 5000) < (seqlen(t)/3)) m = -1; }}}
if (m < 0)
 fprintf(sw->f,"%sNot annotated\n",margin);
else
 { a = list[m];
   fprintf(sw->f,"%sMatch with annotated %s %c(%ld,%ld)",
           margin,d->gene[a].species,comp[d->gene[a].comp],
           d->gene[a].start,d->gene[a].stop);
   w = gene_mismatch_report(d,d->gene+a,t,mreport,sw);
   if (w > 0) fprintf(sw->f," * %s",mreport);
   fputc('\n',sw->f); }
while (++m < n)
 { a = list[m];
   fprintf(sw->f,"%sOverlap with annotated %s %c(%ld,%ld)\n",
           margin,d->gene[a].species,comp[d->gene[a].comp],
           d->gene[a].start,d->gene[a].stop); }
fputc('\n',sw->f);
}

void disp_gene_set(data_set *d, int nt, csw *sw)
{ int i,j,n,vsort[NT],*sort;
  char m[MATX][MATY],s[20];
  gene *t;
  FILE *f = sw->f;
  if (nt <= NT)
   sort = vsort;
  else
   { sort = (int *)malloc(nt*sizeof(int));
     if (sort == NULL)
      { fprintf(stderr,"Not enough memory to sort genes\n");
        exit(1); }}
  n = gene_sort(d,nt,sort,sw);
  j = sw->tmrna_struct[54];
  for (i = 55; i <= 60; i++) j += sw->tmrna_struct[i];
  if (j != ((sw->tmrna_struct[0] << 4) + 9)) return;
  if (sw->libflag < 2)
   { if (n > 0)
      for (j = 0; j < n;)
       { i = sort[j++];
         t = ts + i;
         t->energy = nenergy(t,sw);
         switch(t->genetype)
          { case tRNA:
                    init_matrix(m);
                    disp_gene(t,m,sw);
                    sprintf(s,"%d.",j);
                    xcopy(m,0,32,s,length(s));
                    disp_matrix(f,m,MATY);
	                if (sw->matchacceptor)
                     if (iamatch(d,t,sw) == 0)
		              { fprintf(f,"    Iso-acceptor mismatch\n");
			            sw->iamismatch++; }
                    if (sw->annotated)
                     annotation_overlap_check(d,t,"    ",sw);
                    overlap(d,sort,n,i,sw);
                    if (sw->secstructdisp & 2) disp_trna_bracket_notation(f,t,sw);
                    if (sw->secstructdisp & 4) disp_gene_SVG(t,m,sw);
                    if (sw->seqdisp) disp_seq(f,t,sw);
                    if (t->nintron > 0) disp_intron(f,t,sw);
                    if (sw->energydisp > 1) trna_score(f,t);
                    break;
            case tmRNA:
                    if ((sw->secstructdisp & 1) || (sw->secstructdisp & 4))
                     { init_matrix(m);
                       disp_gene(t,m,sw); }
                    if (sw->secstructdisp & 1)
                     { sprintf(s,"%d.",j);
                       xcopy(m,0,32,s,length(s));
                       disp_matrix(f,m,MATY);
                       if (sw->annotated)
                        annotation_overlap_check(d,t,"    ",sw); }
                    else
                     { fprintf(f,"\n%d.\n",j);
                       disp_location(t,sw,"Location");
                       if (sw->reportpseudogenes)
                        if (pseudogene(t,sw))
                         fprintf(f,"Possible Pseudogene\n");
                       if (sw->energydisp)
                        fprintf(f,"Score = %g\n",t->energy);
                       if (sw->annotated)
                        annotation_overlap_check(d,t,"",sw); }
                    overlap(d,sort,n,i,sw);
                    if (t->asst == 0) disp_tmrna_seq(f,t,sw);
                    else disp_tmrna_perm_seq(f,t,sw);
                    if (sw->secstructdisp & 4) disp_gene_SVG(t,m,sw);
                    if (sw->energydisp > 1) tmrna_score(f,t,sw);
                    break;
            case CDS:
                    fprintf(f,"\n%d.\nCDS gene\n",j);
                    disp_location(t,sw,"Location");
                    if (sw->annotated)
                     annotation_overlap_check(d,t,"",sw);
                    overlap(d,sort,n,i,sw);
                    disp_cds(f,t,sw);
                    break;
                  }
         if (sw->libflag > 0) write_to_library(f,t,sw); }
     else
      if (*(d->seqname) != '\0')
       fprintf(f,"\nNothing found in %s\n\n\n",d->seqname);
      else
       fprintf(f,"\nNothing found\n\n\n"); }
  else
   { if (n > 0)
      for (i = 0; i < n; i++)
       write_to_library(f,ts + sort[i],sw); }
  disp_energy_stats(d,nt,sw);
  if (d->datatype == GENBANK) disp_match(d,sort,n,sw);
  if (nt > NT) free((void *)sort); }


void batch_gene_set(data_set *d, int nt, csw *sw)
{ int i,j,n,vsort[NT],nspaces,caps,*sort;
  gene *t;
  FILE *f = sw->f;
  if (nt <= NT)
   sort = vsort;
  else
  {  sort = (int *)malloc(nt*sizeof(int));
     if (sort == NULL)
      { fprintf(stderr,"Not enough memory to sort genes\n");
        exit(1); }}
  n = gene_sort(d,nt,sort,sw);
  j = sw->tmrna_struct[54];
  for (i = 55; i <= 60; i++) j += sw->tmrna_struct[i];
  if (j != ((sw->tmrna_struct[0] << 4) + 9)) return;
  if (sw->libflag < 2)
   if (sw->batch >= 2)
    { nspaces = (sw->batch & 0x4);
      caps = (sw->batch & 0x10);
      if (sw->batch & 0x8) 
       for (i = 0; i < n; i++)
        disp_fasta_seq(f,ts + sort[i],d->ns+1,i+1,nspaces,caps,sw);
      else
       for (i = 0; i < n; i++)
        disp_fasta_seq(f,ts + sort[i],0,0,nspaces,caps,sw); }
   else
    { if (n == 1) fprintf(f,"1 gene found\n");
      else fprintf(f,"%d genes found\n",n);
      for (j = 0; j < n; j++)
       { fprintf(f,"%-3d ",j+1);
         t = ts + sort[j];
         t->energy = nenergy(t,sw);
         switch(t->genetype)
          { case tRNA:  disp_batch_trna(f,t,sw);
                        break;
            case tmRNA: disp_batch_tmrna(f,t,sw);
                        break;
            case srpRNA:disp_batch_srprna(f,t,sw);
                        break;
            case CDS:   disp_batch_cds(f,t,sw);
                        break;
            default:    break; }}}
  if (sw->libflag > 0)
   { for (i = 0; i < n; i++)
      write_to_library(f,ts + sort[i],sw); }
  batch_energy_stats(d,nt,sw);
  if (nt > NT) free((void *)sort); }


void remove_overlapping_trna(data_set *d, int nt, csw *sw)
{ int i,n,ioverlay;
  long a,b,c,e,len,leni,overlap,psmax;
  char s1[80],s2[80];
  gene *t,*ti;
  static long proximity = 7*MINCTRNALEN/10;
  psmax = d->psmax;
  ioverlay = sw->ioverlay;
  for (n = 0; n < nt; n++)
   { t = ts + n;
     if (t->genetype != tRNA) continue;
     if (t->energy < 0.0) continue;
     if (t->nintron <= 0) continue;
     a = t->start;
     b = t->stop;
     if (b < a) b += psmax;
     len = b - a;
     for (i = 0; i < nt; i++)
      { if (i == n) continue;
        ti = ts + i;
        if (ti->genetype != tRNA) continue;
        if (ti->comp != t->comp) continue;
        if (ti->energy < 0.0) continue;
        c = ti->start;
        e = ti->stop;
        if (e < c) e += psmax;
        leni = e - c;
        if (ioverlay)
         { if ((2*len) > (5*leni)) continue;
           if ((2*leni) > (5*len)) continue; }
        overlap = (a >= c)?((b >= e)?e-a:len):((b >= e)?len:b-c);
        if (overlap >= proximity)
         if (t->energy < ti->energy)
          { if (sw->verbose)
             { fprintf(stderr,"Removing %s at %s",name(t,s1,0,sw),position(s2,t,sw));
               if (sw->energydisp) fprintf(stderr," (%g)",nenergy(t,sw));
               fprintf(stderr,"\n"); }
             t->energy = -1.0;
             break; }}}
  for (n = 0; n < (nt-1); n++)
   { t = ts + n;
     if (t->genetype != tRNA) continue;
     if (t->energy < 0.0) continue;
     a = t->start;
     b = t->stop;
     if (b < a) b += psmax;
     len = b - a;
     for (i = n + 1; i < nt; i++)
      { ti = ts + i;
        if (ti->genetype != tRNA) continue;
        if (ti->comp != t->comp) continue;
        if (ti->energy < 0.0) continue;
        c = ti->start;
        e = ti->stop;
        if (e < c) e += psmax;
        leni = e - c;
        if (ioverlay)
         { if ((2*len) > (5*leni)) continue;
           if ((2*leni) > (5*len)) continue; }
        overlap = (a >= c)?((b >= e)?e-a:len):((b >= e)?len:b-c);
        if (overlap >= proximity)
         if (t->energy < ti->energy)
          { if (sw->verbose)
             { fprintf(stderr,"Removing %s at %s",name(t,s1,0,sw),position(s2,t,sw));
               if (sw->energydisp) fprintf(stderr," (%g)",nenergy(t,sw));
               fprintf(stderr,"\n"); }
            t->energy = -1.0;
            break; }
        else if (ti->energy < t->energy)
          { if (sw->verbose)
             { fprintf(stderr,"Removing %s at %s",name(ti,s1,0,sw),position(s2,ti,sw));
               if (sw->energydisp) fprintf(stderr," (%g)",nenergy(ti,sw));
               fprintf(stderr,"\n"); }
            ti->energy = -1.0; }}}}



void iopt_fastafile(data_set *d, csw *sw)
{ int i,nt,flag,len,aragorn,anticodon;
  int *s,*sf,*se,*sc,*swrap;
  int seq[2*LSEQ+WRAP+1],cseq[2*LSEQ+WRAP+1],wseq[2*WRAP+1];
  long gap,start,rewind,drewind,psmax,tmaxlen,vstart,vstop;
  double sens,sel1,sel2;
  char c1,c2,c3;
  static char trnatypename[3][25] =
   { "Metazoan mitochondrial","Cytosolic","Mammalian mitochondrial" };
  static char genecodename[NGENECODE][50] =
   { "composite Metazoan mitochondrial",
     "standard",
     "vertebrate mitochondrial",
     "yeast mitochondrial",
     "mold/protozoan/Coelenterate mitochondrial",
     "invertebrate mitochondrial",
     "Ciliate",
     "deleted -> standard",
     "deleted -> standard",
     "Echinoderm/Flatworm mitochondrial",
     "Euplotid",
     "bacterial/plant chloroplast",
     "alternative yeast",
     "Ascidian mitochondrial",
     "alternative flatworm mitochondrial",
     "Blepharisma",
     "Chlorophycean mitochondrial",
     "deleted -> standard",
     "deleted -> standard",
     "deleted -> standard",
     "deleted -> standard",
     "Trematode mitochondrial",
     "Scenedesmus obliquus mitochondrial",
     "Thraustochytrium mitochondrial",
     "Pterobranchia mitochondrial",
     "Gracilibacteria",
     "Pachysolen tannophilus",
     "Karyorelict",
     "Condylostoma",
     "Mesodinium",
     "Peritrich",
     "Blastocrithidia",
     "vacant -> standard",
     "Cephalodiscidae mitochondrial UAA-Tyr"
   };
  FILE *f = sw->f;
  init_tmrna(f,sw);
  aragorn = (sw->trna || sw->tmrna || sw->cds || sw->srprna);
  fprintf(f,"\nPlease reference the following paper");
  if (aragorn && sw->mtrna) fputc('s',f);
  fprintf(f," if you use this\n");
  fprintf(f,"program as part of any published research.\n\n");
  if (aragorn)
   { fprintf(f,"Laslett, D. and Canback, B. (2004) ARAGORN, a\n");
     fprintf(f,"program for the detection of transfer RNA and\n");
     fprintf(f,"transfer-messenger RNA genes in nucleotide sequences.\n");
     fprintf(f,"Nucleic Acids Research, 32;11-16.\n\n"); }
  if (sw->mtrna)
   { fprintf(f,"Laslett, D. and Canback, B. (2008) ARWEN: a\n");
     fprintf(f,"program to detect tRNA genes in metazoan mitochondrial\n");
     fprintf(f,"nucleotide sequences\n");
     fprintf(f,"Bioinformatics, 24(2); 172-175.\n\n\n"); }
  fputc('\n',f);
  if (sw->mtrna)
   { fprintf(f,"Searching for %s tRNA genes\n",trnatypename[sw->discrim]);
     if (!sw->tvloop)
      fprintf(f,"TV replacement loop tRNA genes not detected\n"); }
  else
   if (sw->trna)
    { fprintf(f,"Searching for tRNA genes");
      if (sw->maxintronlen > 0) fprintf(f," with introns in anticodon loop");
      else fprintf(f," with no introns");
      fputc('\n',f);
      if (sw->maxintronlen > 0)
       { fprintf(f,"Intron length from %d to %d bases\n",
               sw->minintronlen,sw->maxintronlen);
         if (sw->ifixedpos)
	      { fprintf(f,"Intron position fixed between positions 37 and 38\n");
            fprintf(f,"on C-loop (one base after anticodon)\n"); }
         if (sw->ioverlay)
          fprintf(f,"Allowing overlay of long tRNA genes\n"); }}
  if (sw->tmrna)
    fprintf(f,"Searching for tmRNA genes\n");
  if (sw->linear)
   fprintf(f,"Assuming linear topology, search will not wrap around ends\n");
  else
   fprintf(f,"Assuming circular topology, search wraps around ends\n");
  if (sw->both == 2)
   fprintf(f,"Searching both strands\n");
  else
   if (sw->both == 1)
    fprintf(f,"Searching complementary (antisense) strand only\n");
   else
    fprintf(f,"Searching single (sense) strand only\n");
  if (sw->mtrna)
   if (sw->mtcompov)
    fprintf(f,"Reporting overlapping candidates on opposite strands\n"); 
  if ((sw->mtrna) || (sw->trna) || (sw->tmrna))
   { fprintf(f,"Using %s genetic code\n",genecodename[sw->geneticcode]);
     if (sw->ngcmod > 0)
      { fprintf(f,"Specified modifications:\n");
        for (i = 0; i < sw->ngcmod; i++)
         { anticodon = sw->gcmod[i];
           c1 = cpbase(Thymine - (anticodon & 0x3));
           c2 = cpbase(Thymine - ((anticodon >> 2) & 0x3));
           c3 = cpbase(Thymine - ((anticodon >> 4) & 0x3));
           fprintf(f,"%c%c%c = %s\n",c1,c2,c3,
                   aaname[aamap[sw->geneticcode][anticodon]]); }}}
  fputc('\n',f);
  fputc('\n',f);
  rewind = MAXTAGDIST + 20;
  if (sw->trna | sw->mtrna)
   { tmaxlen = MAXTRNALEN + sw->maxintronlen;
     if (rewind < tmaxlen) rewind = tmaxlen; }
  if (sw->tmrna)
   if (rewind < MAXTMRNALEN) rewind = MAXTMRNALEN;
  if (sw->peptide)
   if (sw->tagthresh >= 5)
    if (rewind < TSWEEP) rewind = TSWEEP;
  sw->loffset = rewind;
  sw->roffset = rewind;
  drewind = 2*rewind;
  d->ns = 0;
  d->nf = 0;
  d->nextseq = 0L;
  d->nextseqoff = 0L;
  while (d->nextseq >= 0L)
   { d->seqstart = d->nextseq;
     d->seqstartoff = d->nextseqoff;
     if (!seq_init(d,sw)) break;
     psmax = d->psmax;
     if (sw->verbose)
      { fprintf(stderr,"%s\n",d->seqname);
        fprintf(stderr,"%ld nucleotides in sequence\n",psmax);
        fprintf(stderr,"Mean G+C content = %2.1f%%\n",100.0*d->gc);
        if ((sw->mtrna) || (sw->trna) || (sw->tmrna))
         { fprintf(stderr,"Using %s genetic code\n",genecodename[sw->geneticcode]);
           if (sw->ngcmod > 0)
            { fprintf(stderr,"Specified modifications:\n");
              for (i = 0; i < sw->ngcmod; i++)
               { anticodon = sw->gcmod[i];
                 c1 = cpbase(Thymine - (anticodon & 0x3));
                 c2 = cpbase(Thymine - ((anticodon >> 2) & 0x3));
                 c3 = cpbase(Thymine - ((anticodon >> 4) & 0x3));
                 fprintf(stderr,"%c%c%c = %s\n",c1,c2,c3,
                         aaname[aamap[sw->geneticcode][anticodon]]); }}}}
     fprintf(f,"%s\n",d->seqname);
     fprintf(f,"%ld nucleotides in sequence\n",psmax);
     fprintf(f,"Mean G+C content = %2.1f%%\n",100.0*d->gc);
     init_gene(0,NT);
     nt = 0;
     flag = 0;
     start = 1L;
     se = seq;
     if (sw->linear)
      { for (i = 0; i < rewind; i++) *se++ = NOBASE;
        start -= rewind; }
     else
      { if (psmax <= drewind)
         { gap = drewind - psmax;
           sc = se + gap;
           while (se < sc) *se++ = NOBASE;
           swrap = wseq;
           sc = se + psmax;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }
           sc = swrap + gap;
           while (swrap < sc) *swrap++ = NOBASE;
                   swrap = wseq;
                   sc = swrap + psmax;
                   while (swrap < sc) *se++ = *swrap++;
                   swrap = wseq;
                   sc = swrap + drewind;
                   while (swrap < sc) *se++ = *swrap++;
                   sw->loffset = drewind;
                   sw->roffset = drewind;
                   start -= drewind;
                   flag = 1;
                   goto SH; }
        else
         { swrap = wseq;
           sc = seq + drewind;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }}}
     sc = seq + LSEQ;
     NX:
     while (se < sc)
      { if (d->ps >= psmax)
         { if (sw->linear)
            for (i = 0; i < rewind; i++) *se++ = NOBASE;
           else
            { sc = wseq + drewind;
              swrap = wseq;
              while (swrap < sc) *se++ = *swrap++; }
           flag = 1;
           break; }
        else *se++ = move_forward(d); }
     SH:
     len = (int)(se - seq);
     if (sw->verbose)
      { vstart = sq(start + sw->loffset);
        vstop = sq(start + len - sw->roffset - 1);
        if (vstop < vstart)
         { fprintf(stderr,"Searching from %ld to %ld\n",vstart,psmax);
           fprintf(stderr,"Searching from 1 to %ld\n",vstop); }
        else
         fprintf(stderr,"Searching from %ld to %ld\n",vstart,vstop); }
     if (sw->both != 1)
      { sw->start = start;
        sw->comp = 0;
        nt = tmioptimise(d,seq,len,nt,sw); }
     if (sw->both > 0)
      { sense_switch(seq,cseq,len);
        sw->start = start+len;
        sw->comp = 1;
        nt = tmioptimise(d,cseq,len,nt,sw); }
     if (!flag)
      { s = seq;
        sf = se - drewind;
        se = seq + drewind;
        while (s < se) *s++ = *sf++;
        start += len - drewind;
        goto NX; }
     if (nt < 1) d->nf++;
     if (sw->maxintronlen > 0) remove_overlapping_trna(d,nt,sw);
     if (sw->updatetmrnatags) update_tmrna_tag_database(ts,nt,sw);
     disp_gene_set(d,nt,sw);
     if (sw->verbose) fprintf(stderr,"%s\nSearch Finished\n\n",d->seqname);
     d->ns++; }
  if (d->ns > 1)
   { fprintf(f,"\n\n%d sequences searched\n",d->ns);
     if (sw->trna | sw->mtrna)
      { fprintf(f,"Total tRNA genes = %d\n",sw->ngene[tRNA]);
	    if (sw->matchacceptor)
         fprintf(f,"Total iso-acceptor mismatches = %d\n",sw->iamismatch); }
     if (sw->tmrna) fprintf(f,"Total tmRNA genes = %d\n",sw->ngene[tmRNA]); 
     if (sw->reportpseudogenes)
      if (sw->nps > 0)
       fprintf(f,"Total number of possible pseudogenes = %d\n",sw->nps);
     if (d->nf > 0)
      { sens = 100.0*(d->ns - d->nf)/d->ns; 
        fprintf(f,"Nothing found in %d sequences (%.2lf%% sensitivity)\n",d->nf,sens); }
     if (sw->annotated)
      { if (sw->trna | sw->mtrna) 
         { fprintf(f,"\nTotal number of annotated tRNA genes = %d\n",
                   sw->nagene[tRNA]);
           fprintf(f,"Total number of annotated tRNA genes not detected = %d\n",sw->nafn[tRNA]);
           fprintf(f,"Total number of unannotated tRNA genes detected = %d\n",sw->nafp[tRNA]);
           fprintf(f,"Total number of unannotated DRL tRNA genes detected = %d\n",
                   sw->natfpd);
           fprintf(f,"Total number of unannotated TVRL tRNA genes detected = %d\n",
                   sw->natfptv);
           fprintf(f,"Total annotated sequence length = %ld bases\n",sw->nabase);
           sens = (sw->nagene[tRNA] > 0)?
                         100.0*(double)(sw->nagene[tRNA] - sw->nafn[tRNA])/
                         (double)sw->nagene[tRNA]:0.0;
           sel1 = (sw->nagene[tRNA] > 0)?
                         100.0*(double)(sw->nafp[tRNA])/
                         (double)sw->nagene[tRNA]:0.0;
           sel2 = (sw->nabase > 0)?
                         1000000.0*(double)(sw->nafp[tRNA])/
                         (double)sw->nabase:0.0;
           fprintf(f,"Sensitivity = %lg%%\n",sens);
           fprintf(f,"Selectivity = %lg%% or %lg per megabase\n\n",sel1,sel2); }
        if (sw->tmrna) 
         { fprintf(f,"\nTotal number of annotated tmRNA genes = %d\n",
                   sw->nagene[tmRNA]);
           fprintf(f,"Total number of annotated tmRNA genes not detected = %d\n",sw->nafn[tmRNA]);
           fprintf(f,"Total number of unannotated tmRNA genes detected = %d\n",sw->nafp[tmRNA]);
           fprintf(f,"Total annotated sequence length = %ld bases\n",sw->nabase);
           sens = (sw->nagene[tmRNA] > 0)?
                         100.0*(double)(sw->nagene[tmRNA] - sw->nafn[tmRNA])/
                         (double)sw->nagene[tmRNA]:0.0;
           sel1 = (sw->nagene[tmRNA] > 0)?
                         100.0*(double)(sw->nafp[tmRNA])/
                         (double)sw->nagene[tmRNA]:0.0;
           sel2 = (sw->nabase > 0)?
                         1000000.0*(double)(sw->nafp[tmRNA])/
                         (double)sw->nabase:0.0;
           fprintf(f,"Sensitivity = %lg%%\n",sens);
           fprintf(f,"Selectivity = %lg%% or %lg per Megabase\n\n",sel1,sel2); }
        if (sw->cds) 
         { fprintf(f,"\nTotal number of annotated CDS genes = %d\n",
                   sw->nagene[CDS]);
           fprintf(f,"Total number of annotated CDS genes not detected = %d\n",sw->nafn[CDS]);
           fprintf(f,"Total number of unannotated CDS genes detected = %d\n",sw->nafp[CDS]);
           fprintf(f,"Total annotated sequence length = %ld bases\n",sw->nabase);
           sens = (sw->nagene[CDS] > 0)?
                         100.0*(double)(sw->nagene[CDS] - sw->nafn[CDS])/
                         (double)sw->nagene[CDS]:0.0;
           sel1 = (sw->nagene[CDS] > 0)?
                         100.0*(double)(sw->nafp[CDS])/
                         (double)sw->nagene[CDS]:0.0;
           sel2 = (sw->nabase > 0)?
                         1000000.0*(double)(sw->nafp[CDS])/
                         (double)sw->nabase:0.0;
           fprintf(f,"Sensitivity = %lg%%\n",sens);
           fprintf(f,"Selectivity = %lg%% or %lg per Megabase\n",sel1,sel2);
           sens = (sw->lacds > 0)?
                         100.0*(double)sw->ldcds/(double)sw->lacds:0.0;
           fprintf(f,"Length sensitivity = %lg%%\n\n",sens); }
      } }
  if (sw->updatetmrnatags) report_new_tmrna_tags(sw);
}


void bopt_fastafile(data_set *d, csw *sw)
{ int i,nt,flag,len;
  int *s,*sf,*se,*sc,*swrap;
  int seq[2*LSEQ+WRAP+1],cseq[2*LSEQ+WRAP+1],wseq[2*WRAP+1];
  long gap,start,rewind,drewind,psmax,tmaxlen,vstart,vstop;
  double sens;
  FILE *f = sw->f;
  rewind = MAXTAGDIST + 20;
  if (sw->trna | sw->mtrna)
   { tmaxlen = MAXTRNALEN + sw->maxintronlen;
     if (rewind < tmaxlen) rewind = tmaxlen; }
  if (sw->tmrna)
   if (rewind < MAXTMRNALEN) rewind = MAXTMRNALEN;
  if (sw->peptide)
   if (sw->tagthresh >= 5)
    if (rewind < TSWEEP) rewind = TSWEEP;
  sw->loffset = rewind;
  sw->roffset = rewind;
  drewind = 2*rewind;
  d->ns = 0;
  d->nf = 0;
  d->nextseq = 0L;
  d->nextseqoff = 0L;
  while (d->nextseq >= 0L)
   { d->seqstart = d->nextseq;
     d->seqstartoff = d->nextseqoff;
     if (!seq_init(d,sw)) break;
     psmax = d->psmax;
     if (sw->verbose)
      { fprintf(stderr,"%s\n",d->seqname);
        fprintf(stderr,"%ld nucleotides in sequence\n",psmax);
        fprintf(stderr,"Mean G+C content = %2.1f%%\n",100.0*d->gc); }
     if (sw->batch < 2) fprintf(f,">%s\n",d->seqname);
     init_gene(0,NT);
     nt = 0;
     flag = 0;
     start = 1L;
     se = seq;
     if (sw->linear)
      { for (i = 0; i < rewind; i++) *se++ = NOBASE;
        start -= rewind; }
     else
      { if (psmax <= drewind)
         { gap = drewind - psmax;
           sc = se + gap;
           while (se < sc) *se++ = NOBASE;
           swrap = wseq;
           sc = se + psmax;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }
           sc = swrap + gap;
           while (swrap < sc) *swrap++ = NOBASE;
                   swrap = wseq;
                   sc = swrap + psmax;
                   while (swrap < sc) *se++ = *swrap++;
                   swrap = wseq;
                   sc = swrap + drewind;
                   while (swrap < sc) *se++ = *swrap++;
                   sw->loffset = drewind;
                   sw->roffset = drewind;
                   start -= drewind;
                   flag = 1;
                   goto SH; }
        else
         { swrap = wseq;
           sc = seq + drewind;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }}}
     sc = seq + LSEQ;
     NX:
     while (se < sc)
      { *se++ = move_forward(d);
        if (d->ps >= psmax)
         { if (sw->linear)
            for (i = 0; i < rewind; i++) *se++ = NOBASE;
           else
            { sc = wseq + drewind;
              swrap = wseq;
              while (swrap < sc) *se++ = *swrap++; }
           flag = 1;
           break; }}
         SH:
     len = (int)(se - seq);
     if (sw->verbose)
      { vstart = sq(start + sw->loffset);
        vstop = sq(start + len - sw->roffset - 1);
        if (vstop < vstart)
         { fprintf(stderr,"Searching from %ld to %ld\n",vstart,psmax);
           fprintf(stderr,"Searching from 1 to %ld\n",vstop); }
        else
         fprintf(stderr,"Searching from %ld to %ld\n",vstart,vstop); }
     if (sw->both != 1)
      { sw->start = start;
        sw->comp = 0;
        nt = tmioptimise(d,seq,len,nt,sw); }
     if (sw->both > 0)
      { sense_switch(seq,cseq,len);
        sw->start = start+len;
        sw->comp = 1;
        nt = tmioptimise(d,cseq,len,nt,sw); }
     if (!flag)
      { s = seq;
        sf = se - drewind;
        se = seq + drewind;
        while (s < se) *s++ = *sf++;
        start += len - drewind;
        goto NX; }
     if (nt < 1) d->nf++;
     if (sw->maxintronlen > 0) remove_overlapping_trna(d,nt,sw);
     if (sw->updatetmrnatags) update_tmrna_tag_database(ts,nt,sw);
     batch_gene_set(d,nt,sw);
     if (sw->verbose) fprintf(stderr,"%s\nSearch Finished\n\n",d->seqname);
     d->ns++; }
  if ((d->ns > 1) && (sw->batch < 2))
   { fprintf(f,">end \t%d sequences",d->ns);
     if (sw->trna || sw->mtrna) fprintf(f," %d tRNA genes",sw->ngene[tRNA]);
     if (sw->tmrna) fprintf(f," %d tmRNA genes",sw->ngene[tmRNA]);
     if (d->nf > 0)
      { sens = 100.0*(d->ns - d->nf)/d->ns; 
        fprintf(f,", nothing found in %d sequences, (%.2lf%% sensitivity)",d->nf,sens); }
     fputc('\n',f); }
  if (sw->updatetmrnatags) report_new_tmrna_tags(sw);
}


void aragorn_help_menu()
{
int h;
for (h = 0; h < NHELPLINE; h++) printf("%s\n",helpmenu[h]);
}

void error_report(int n, char *s)
{ switch(n)
   { case 0: fprintf(stderr,
              "-%s not recognised, type aragorn -h for help\n",s);
             break;
     case 1: fprintf(stderr,
              "-%s not understood, type aragorn -h for help\n",s);
             break;
     case 2: fprintf(stderr,"Could not open %s\n",s);
             break;
     case 3: fprintf(stderr,
             "No sequence file specified, type aragorn -h for help\n");
             break;
     case 4: fprintf(stderr,"Don't know genetic code %s\n",s);
             break;
     case 5: fprintf(stderr,"Too many genetic code modifications (max=%d)\n",
                     MAXGCMOD);
             break;
     default: break; }
  exit(0); }


void process_genecode_switch(char *s, csw *sw)
{ int i,m,lmax,len[NGENECODE],anticodon,b[3];
  long l;
  char c,*ss,*se;
  static char genecodetag[NGENECODE][10] =
   { "MET",
     "STD","VERT","YEAST","PROT","INVERT",
     "CILIATE","DELETED","DELETED","FLATWORM","EUPLOT",
     "BACT","ALTYEAST","ASCID","ALTFLAT","BLEP",
     "CHLOROPH","DELETED","DELETED","DELETED","DELETED",
     "TREM","SCEN","THRAUST","PTERO","GRAC",
     "PACH","KARY","COND","MESO","PERI","BLAST","VACANT","CEPH" };
  sw->geneticcode = STANDARD;
  sw->gcfix = 1;
  c = *s;
  if (c >= '0')
   if (c <= '9')
    { lconvert(s,&l);
      i = (int)l;
      if ((i >= 0) && (i < NGENECODE)) sw->geneticcode = i;
      goto MOD; }
  for (i = 0; i < NGENECODE; i++) 
   { len[i] = 0;
     ss = s;
     se = genecodetag[i];
     while (c == *ss++)
      { if (upcasec(c) != *se++) break;
        len[i]++; }}
  m = -1;
  lmax = 0;
  i = -1;
  while (++i < NGENECODE)
   if (len[i] > lmax)
    { m = i;
      lmax = len[i]; }
  if (m >= 0) sw->geneticcode = m;
  else error_report(4,s);
  MOD:
  sw->ngcmod = 0;
  ss = s;
  while (ss = strpos(ss,","))
   { if (sw->ngcmod >= MAXGCMOD) error_report(5,NULL);
     ss++;
     for (i = 0; i < 3; i++)
      { b[i] = Adenine;
        c = upcasec(ss[i]);
        if (c == 'C') b[i] = Cytosine;
        if (c == 'G') b[i] = Guanine;
        if (c == 'T') b[i] = Thymine;
        if (c == 'U') b[i] = Thymine; }
      anticodon = ((Thymine - b[2])<<4) + ((Thymine - b[1])<<2) + (Thymine - b[0]);
     if (!(se = strpos(ss,"="))) break;
     se++;
     for (i = 0; i < NAMINOACID; i++)
      if (upcasec(se[0]) == upcasec(aaname[i][0]))
       if (upcasec(se[1]) == upcasec(aaname[i][1]))
        if (upcasec(se[2]) == upcasec(aaname[i][2]))
         { aamap[sw->geneticcode][anticodon] = i;
           sw->gcmod[sw->ngcmod] = anticodon;
           sw->ngcmod++;
           break; } }}
    


void change_thresholds(csw *sw, double psthresh)
{
sw->threshlevel = psthresh;
sw->cdsthresh *= psthresh;
sw->srpthresh *= psthresh;
sw->tmrnathresh *= psthresh;
sw->mtdtthresh *= psthresh;
sw->mttthresh *= psthresh;
sw->mtdthresh *= psthresh;
sw->trnathresh *= psthresh;
}


int main(int z, char *v[])
{ int i,lv,filecounter;
  long l;
  double psthresh;
  char c1,c2,c3,c4,*s;
  data_set d;
  static csw sw =
   { {"tRNA","tmRNA","","","CDS","overall"},
     NULL,0,0,0,0,0,0,0,0,1,0,0,
     STANDARD,0,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},0,METAZOAN_MT,
     1,0,5,5,1,0,0,0,2,0,0,0,0,0,0,3,0,2,1,1,0,0,0,0,0,0,0,0,1,
     0,0,0,0,0,0,0,{0,0,0,0,0,0},0,0,0,0,NTAG,10,30,
     {0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},0,0,0,0,0L,
     100.0,1.0,tRNAthresh,4.0,29.0,26.0,7.5,8.0,
     mtRNAtthresh,mtRNAdthresh,mtRNAdtthresh,-7.9,-6.0,
     tmRNAthresh,14.0,10.0,25.0,9.0,srpRNAthresh,CDSthresh,
     {tRNAthresh,tmRNAthresh,srpRNAthresh,0.0,CDSthresh },
     { 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       10, 65, 82, 65, 71, 79, 82, 78, 32,
       118, 49, 46, 50, 46, 52, 48, 32, 32, 32,
       68, 101, 97,110, 32, 76, 97, 115, 108,
       101, 116, 116, 10,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       10, TERM }};
  sw.f = stdout;
  d.bugmode = 0;
  filecounter = 0;
  i = 0;
  while (++i < z)
   if (*(v[i]) == '-')
    { lv = length(v[i]);
      if (lv < 2) continue;
      s = v[i] + 1;
      c1 = upcasec(*s);
      c2 = (lv > 2)?upcasec(s[1]):' ';
      c3 = (lv > 3)?upcasec(s[2]):' ';
      c4 = (lv > 4)?upcasec(s[3]):' ';
      switch(c1)
       { case  'E': sw.energydisp = (c2 == 'S')?2:1;
                    break;
         case  'A': if (c2 == '7') sw.extastem = 0;
                    else
                     if (c2 == 'A') sw.matchacceptor = 1;
		             else 
                      if (c2 == 'M')
                       { l = 1L;
                         if (c3 == 'T')
                          { if (lv > 4) 
                             { s = lconvert(s+3,&l);
                               if (l < 1L) l = 1L;
                               sw.trnalenmisthresh = (int)l; }
                            else sw.trnalenmisthresh = 1; }
                         else if (c3 == 'M')
                          { if (lv > 4) 
                             { s = lconvert(s+3,&l);
                               if (l < 1L) l = 1L;
                               sw.tmrnalenmisthresh = (int)l; }
                            else sw.tmrnalenmisthresh = 1; }
                         else if (lv > 3)
                          { s = lconvert(s+2,&l);
                            if (l < 1L) l = 1L;
                            sw.trnalenmisthresh = (int)l;
                            sw.tmrnalenmisthresh = (int)l; }
                         else
                          { sw.trnalenmisthresh = 1;
                            sw.tmrnalenmisthresh = 1; }}
                      else sw.secstructdisp |= 1;
                    break;
         case  'B': if (c2 == 'R') sw.secstructdisp |= 2;
                    else sw.libflag = 1;
                    break;
         case  'X': sw.libflag = 2;
                    break;
         case  'W': if (c2 == 'U')
                     if (c3 == 'N')
                      if (c4 == 'I')
                       { d.bugmode = 1;
                         break; }
                    if (sw.batch < 1) sw.batch = 1;
                    if (c2 == 'A') sw.batchfullspecies = 1;
                    break;
         case  'V': sw.verbose = 1;
                    break;
         case  'S': if (c2 == 'S')
                     { sw.sp1max = 2;
                       sw.sp2min = 1;
                       sw.sp2max = 1;
                       break; }
                    if (c2 == 'E')
                     { if (sw.seqdisp < 1) sw.seqdisp = 1;
                       break; }
                    if ((c2 == 'C') || (c2 == '-'))
                     { sw.both = 1;
                       break; }
                    if ((c2 == 'V') && (c3 == 'G'))
                     { sw.secstructdisp |= 4;
                       break; }
                    sw.both = 0;
                    break;
         case  'F': if (softstrpos(s,"O")) 
                     { sw.batch = 2;
                       if (softstrpos(s,"S")) sw.batch |= 0x4;
                       if (softstrpos(s,"N")) sw.batch |= 0x8;
                       if (softstrpos(s,"C")) sw.batch |= 0x10; }
                    else 
                     { if (softstrpos(s,"C")) sw.seqdisp = 4;
                       else sw.seqdisp = 3; }
                    break;
         case  'D': sw.both = 2;
                    break;
         case  'L': sw.linear = 1;
                    break;
         case  'C': if (c2 == '7') sw.cloop7 = 1;
                    else sw.linear = 0;
                    break;
         case  'J': if (lv > 2)
                     { if (c2 == 'R') sw.aataildiv = 1;
                       if (c3 == '4') sw.aataildisp = 1; }
                    else sw.aataildisp = 1;
                    break;
         case  '1': sw.minintronlen = 10;
                    break;
         case  'I': if (c2 == 'O') { sw.ioverlay = 1; s++; lv--; }
                    else if (c2 == 'F') { sw.ifixedpos = 1; s++; lv--; }
                    else if (c2 == 'R') { sw.ireportminintronlen = 1; s++; lv--; }
                    if (c3 == 'O') { sw.ioverlay = 1; s++; lv--; }
                    else if (c3 == 'F') { sw.ifixedpos = 1; s++; lv--; }
                    else if (c3 == 'R') { sw.ireportminintronlen = 1; s++; lv--; }
                    if (c4 == 'O') { sw.ioverlay = 1; s++; lv--; }
                    else if (c4 == 'F') { sw.ifixedpos = 1; s++; lv--; }
                    else if (c4 == 'R') { sw.ireportminintronlen = 1; s++; lv--; }
                    if (lv > 2) s = lconvert(s+1,&l);
		            else goto IMAX;
                    if (*s == ',')
                     { if (sw.ireportminintronlen == 1)
                        sw.minintronlenreport = (int)l;
                       else
                        sw.minintronlen = (int)l;
                       lconvert(s+1,&l);
                       sw.maxintronlen = (int)l; }
                    else sw.maxintronlen = (int)l;
                    if (sw.maxintronlen > (LSEQ - MAXTRNALEN))
                     sw.maxintronlen = (LSEQ - MAXTRNALEN);
                    if (sw.maxintronlen > MAXINTRONLEN)
                     sw.maxintronlen = MAXINTRONLEN;
                    if ((sw.minintronlen < 0) ||
                     (sw.maxintronlen < sw.minintronlen))
                      error_report(1,v[i]);
                    if ((sw.minintronlenreport < 0) ||
                     (sw.maxintronlen < sw.minintronlenreport))
                      error_report(1,v[i]);
                    break;
                    IMAX:
                    sw.maxintronlen = MAXINTRONLEN;
                    break;
         case  'T': if (c2 == 'V')
                     { sw.tvloop = 0;
                       break; }
                    sw.trna = 1;
                    if (lv > 2)
                     { s = dconvert(s+1,&sw.trnathresh);
                       if (*s == ',') dconvert(s+1,&sw.ttarmthresh); }
                    break;
         case  'M': if (c2 == 'T')
                     { sw.mtrna = 1;
                       if (!sw.gcfix) sw.geneticcode = METAZOAN_MT;
                       if (lv > 3)
                        { s += 2;
                          c3 = upcasec(*s);
                          if (c3 == 'M')
                           { do c3 = upcasec(*++s);
                             while ((c3 == 'A') || (c3 == 'M')
                                     || (c3 == 'L'));
                             sw.tvloop = 0;
                             sw.geneticcode = VERTEBRATE_MT;
                             sw.discrim = MAMMAL_MT; }
                          MTNXTC:
                          if (c3 == 'X')
                           { c3 = upcasec(*++s);
                             sw.mtxdetect = 0;
                             goto MTNXTC; }
                          if (c3 == 'C')
                           { c3 = upcasec(*++s);
                             sw.mtcdsscan = 0;
                             goto MTNXTC; }
                          if (c3 == 'D')
                           { c3 = upcasec(*++s);
                             sw.mtcompov = 1;
                             goto MTNXTC; }
                          if (c3 != '-')
                           if (c3 != '.')
                            if ((c3 < '0') || (c3 > '9'))
                             break;
                          s = dconvert(s,&sw.mtdtthresh);
                          if (*s == ',') s = dconvert(s+1,&sw.mttthresh);
                          if (*s == ',') s = dconvert(s+1,&sw.mtdthresh);
                          if (*s == ',') s = dconvert(s+1,&sw.mttarmthresh);
                          if (*s == ',') dconvert(s+1,&sw.mtdarmthresh); }}
                    else
                     { sw.tmrna = 1;
                       if (c2 == 'U')
                        if (c3 == 'T')
                         { sw.updatetmrnatags = 1;
                           lv -= 2;
                           s += 2; }
                       if (lv > 2) dconvert(s+1,&sw.tmrnathresh); }
                    break;
         case  'P': if (c2 == 'S')
                     { if (c3 != '-')
                        if (c3 != '.')
                         if ((c3 < '0') || (c3 > '9'))
                          { change_thresholds(&sw,PSEUDOGENElevel);
                            break; }
                       psthresh = 1.0;
                       dconvert(s+2,&psthresh);
                       change_thresholds(&sw,psthresh);
                       break; }
                    break;
         case  'G': if (c2 != 'C') break;
                    process_genecode_switch(s+2,&sw);
                    break;
         case  'R': if (c2 == 'N') sw.repeatsn = 1;
                    else 
                     if (c2 == 'P') 
                      { sw.reportpseudogenes = 1;
                        if (lv > 3)
                         dconvert(s+2,&sw.reportpsthresh); }
                     else sw.tmstrict = 0;
                    break;
         case  'Q': sw.showconfig = 0;
                    break;
         case  'H': aragorn_help_menu();
                    exit(0);
         case  'O': if (lv > 2)
                     s++;
                    else
                     { if (++i >= z) break;
                       s = v[i]; }
                    sw.f = fopen(s,"w");
                    if (!sw.f) error_report(2,s);
                    break;
         default:   error_report(0,s); }}
   else
    if (filecounter < 1)
     { d.f = fopen(v[i],"r");
       if (d.f)
        filecounter++;
       else
        error_report(2,v[i]); }
    else
     if (filecounter < 2)
      { sw.f = fopen(v[i],"w");
        if (!sw.f) error_report(2,v[i]);
        filecounter++; }
     else
     error_report(0,v[i]);
  if (filecounter < 1)
   error_report(3,NULL);
  if ((!sw.trna) & (!sw.tmrna))
   { sw.trna = 1;
     sw.tmrna = 1; }
  if (sw.mtrna) sw.trna = 0;
  ts = (gene *)malloc(NT*sizeof(gene));
  if (ts == NULL)
   { fprintf(stderr,"Not enough memory available to store detected genes\n");
     exit(1); }
  sw.genespace = NT;
  if (sw.libflag) fprintf(sw.f,"Library\n");
  if (sw.batch) bopt_fastafile(&d,&sw);
  else iopt_fastafile(&d,&sw); 
  free((void *)ts);
  fclose(d.f);
  if (!sw.batch && sw.showconfig)
   { fprintf(sw.f,"Configuration: ");
     i = -1;
     while (++i < z) fprintf(sw.f,"%s ",v[i]);
     fputc('\n',sw.f); }
  if (sw.f != stdout) fclose(sw.f);
  return(0); }


