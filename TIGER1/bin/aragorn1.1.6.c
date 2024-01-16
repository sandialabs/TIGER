
/*
---------------------------------------------------------------
ARAGORN v1.1 Dean Laslett 
---------------------------------------------------------------

    ARAGORN.C 
    Detects tRNA and tmRNA genes in nucleotide sequences
    Copyright (C) 2003 Dean Laslett 
  
    Please, report bugs and suggestions of improvements to the authors
  
    E-mail: Björn Canbäck: bcanback@thep.lu.se
            Dean Laslett:  1942651@student.murdoch.edu.au 
  
    Version 1.1  April 28th, 2004. 
  
    Please reference the following paper if you use this
    program as part of any published research. 
  
    Laslett, D. and Canback, B. (2004)
    ARAGORN, a program for the detection of transfer RNA and 
    transfer-messenger RNA genes in nucleotide sequences. 
    Nucleic Acids Research, 32;11-16.
  
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
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

		     END OF TERMS AND CONDITIONS

*/




/*
---------------------------------------------------------------
ARAGORN v1.1 Dean Laslett 
---------------------------------------------------------------


aragorn detects tRNA and tmRNA genes.

Usage:
aragorn -v -s -d -c -l -i<min>,<max> -t -m -mt -o
<outfile> <filename>

<filename> is assumed to contain one or more sequences
in FASTA format. Results of the search are printed to
STDOUT. All switches are optional and
case-insensitive.
Unless -i is specified, tRNA genes containing introns
are not detected.

    -m            Search for tmRNA genes only.
    -t            Search for tRNA genes only.
                  By default, both are detected.
    -i            Search for tRNA genes with introns in
                  anticodon loop with maximum length 3000
                  bases. Minimum intron length is 0 bases.
                  Ignored if -m is specified.
    -i<max>       Search for tRNA genes with introns in
                  anticodon loop with maximum length <max>
                  bases. Minimum intron length is 0 bases.
                  Ignored if -m is specified.
    -i<min>,<max> Search for tRNA genes with introns in
                  anticodon loop with maximum length <max>
                  bases, and minimum length <min> bases.
                  Ignored if -m is specified.
    -c            Assume that each sequence has a circular
                  topology. Search wraps around each end.
                  Default setting.
    -l            Assume that each sequence has a linear
                  topology. Search does not wrap.
    -d            Double. Search both strands of each
                  sequence. Default setting.
    -s            Single. Do not search the complementary
                  strand of each sequence.
    -v            Verbose. Prints out information during
                  search to STDERR.
    -a            Print out tRNA domain for tmRNA genes
    -q            Dont print configuration line (which switchs
                  and files were used).
    -w            Print out in Batch mode.
    -O <outfile>  Print output to <outfile>. If <outfile>
                  already exists, it is overwritten.  By default
                  all output goes to stdout.

So far ARAGORN has been compiled and run
successfully on the
following UNIX machines, using the following commands:

SUN                  cc -O3 -o aragorn aragorn.c
Silicon Graphics     cc -O -o aragorn aragorn.c
Linux                gcc -O3 -ffast-math -finline-functions 
                         -o aragorn aragorn.c


*/


#include <stdio.h>
#include <stdlib.h>

#define DLIM            '\n'
#define SEEK_SET        0
#define SEEK_CUR        1
#define SEEK_END        2
#define STRLEN          4001
#define STRLENM1        4000
#define KEYLEN          15
#define INACTIVE        2.0e+35
#define IINACTIVE       2000000001L
#define ITHRESHOLD      2000000000L
#define space(c)        (c==' ')||(c=='\n')||(c=='\t')
#define sq(pos)         ((pos + d->psmax - 1L) % d->psmax) + 1L 
#define itmparam(x,y)   fputc(x,y)

#define TERM            -1
#define NOBASE          4
#define INSERT          5



#define RIGHT   0
#define UP      1
#define LEFT    2
#define DOWN    3
#define UPRIGHT 4

#define MATX 40
#define MATY 34


#define ASTEM1_EXT       0
#define ASTEM2_EXT       9
#define ASTEM2_EXTD      4                   /* <= ASTEM2_EXT */
#define ASTEM2_EXTE      5                   /* ASTEM2_EXT - ASTEM2_EXTD */
#define MINTSTEM_DIST      (17 + ASTEM2_EXT)
#define MAXTSTEM_DIST      (26 + ASTEM2_EXT)
#define MAXDSTEM_DIST      (9 + ASTEM1_EXT)
#define MINDSTEM_DIST      (8 + ASTEM1_EXT)
#define MININTRONLEN    0
#define MAXINTRONLEN    3000
#define MINCTRNALEN     62
#define MAXCTRNALEN     110 
#define MINTRNALEN      (MINCTRNALEN + ASTEM1_EXT + ASTEM2_EXT)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM1_EXT + ASTEM2_EXT)
#define MAXETRNALEN     (MAXTRNALEN + MAXINTRONLEN)
#define VARMAX          25
#define VARMIN          3
#define VARDIFF         22                   /* VARMAX - VARMIN */
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

#define NA   MAXINTRONLEN
#define ND   100
/* #define NT   50 */
#define NT   150
#define NH   2000
#define NTH  3000
#define NC   4000
#define NS   5
#define NTAG 286
#define LSEQ 20000      
#define ATBOND 2.5


typedef struct { char filename[80];
                 FILE *f;
                 char seqname[STRLEN];
                 double gc;
                 long ps;
                 long psmax;
                 long seqstart;
                 long nextseq; } data_set;
                                              

typedef struct { char name[80];
                 int seq[MAXTRNALEN+1];
                 int eseq[MAXETRNALEN+1];
                 int *ps;
                 int nbase;
		 int comp;
                 long start;
		 long stop;
                 int astem1;
		 int astem2;
                 int spacer;
		 int dstem;
		 int dloop;
		 int cstem;
		 int cloop;
		 int intron;
		 int nintron;
                 int anticodon;
                 int var;
		 int tstem;
		 int tloop;
		 int status;
                 double energy;
                 int asst;
                 int tps;
                 int tpe;   } gene;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 double energy; } trna_loop;


typedef struct { int *pos;
                 int *end;
                 int stem;
                 int loop;
                 double energy; } trna_dloop;


typedef struct { int *pos1;
                 int *pos2;
                 double energy; } astem;


typedef struct { char arg[50];
                 FILE *f;
		 int batch;
                 int trna;
                 int tmrna;
		 int mt;
                 int peptide;
                 int tarm;
                 int tagthresh;
                 int tarmlength;
		 int showconfig;
                 int libflag;
                 int verbose;
                 int linear;
                 int both;
                 int energydisp;
                 int trnadomaindisp;
                 int seqdisp;
                 double trnathresh;
                 double tmrnathresh;
                 double tarmthresh;
                 double tcthresh;
                 double tthresh;
                 double dthresh;
                 double tathresh;
                 double tmathresh;
                 double tmcthresh;
                 double tmcathresh;
                 int tmstrict;
                 unsigned int tmrthresh;
                 unsigned int itathresh;
                 int maxintronlen;
                 int minintronlen;
		 double mttthresh;
		 double mtdthresh;
		 double mtdtthresh;
		 double mttarmthresh;
		 double mtdarmthresh;
                 int ngene[NS];
                 int loffset;
                 int roffset;
                 long start;
                 int comp;
                 int space;
                 int tmrna_struct[200];
               } csw;


/* TOOLS */

char upcasec(char c)
{ return((c >= 'a')?c-32:c); }


int length(char *s)
{ int i=0;
  while (*s++) i++;
  return(i); }


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
  *r = result;
  return(s); }

char *iconvert(char *s, int *r)
{ static char zero='0',nine='9';
  int sgn;
  int result;
  char c;
  sgn = 1;
  result = 0;
  if ((c = *s) == '-')
   { sgn = -1;
     c = *++s; }
  else if (c == '+') c= *++s;
  if (c >= zero)
   if (c <= nine)
    { result = (int)(c - zero);
      while ((c = *++s) >= zero)
       { if (c > nine) break;
         result = result*10 + (int)(c - zero); }}
  *r = result * sgn;
  return(s); }



char *copy(char *from, char *to)
{ while (*to++ = *from++);
  return(--to);  }


char *copy3cr(char *from, char *to, int n)
{ while (*to = *from++)
   { if (*to == DLIM)
      { *to = '\0';
        break; }
     if (--n <= 0) break;
     to++; }
  return(to); }

/* LIBRARY */

long get_seqname(data_set *d)
{ int ic;
  char c,*s;
  FILE *f;
  f = d->f;
  fseek(f,d->seqstart,SEEK_SET);
  do { if ((ic = getc(f)) == EOF) return(-1L);
       c = (char)ic; }
  while (space(c));
  if (c != '>')
   { *d->seqname = '\0';
     fseek(f,d->seqstart,SEEK_SET);
     return(d->seqstart); }
  else
   { if (!fgets(d->seqname,STRLENM1,f)) return(-1L);
     s = d->seqname + length(d->seqname);
     do { if (s < d->seqname) break;
          c = *--s; }
     while ((c == '\n') || (c == '\r'));
     *++s = '\0';
     return(ftell(f)); }}


int move_forward(data_set *d)
{ int ic;
  static int map[256] =
  { -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,NOBASE,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-2,-4,-4, 0,NOBASE, 1,NOBASE,-4,-4, 2,NOBASE,-4,-4,NOBASE,
    -4,NOBASE,NOBASE,-4,-4,-4,NOBASE,NOBASE, 3, 3,NOBASE,NOBASE,-4,
    NOBASE,-4,-4,-4,-4,INSERT,NOBASE,-4, 0,NOBASE, 1,NOBASE,-4,-4, 2,
    NOBASE,-4,-4,NOBASE,-4,NOBASE,NOBASE,-4,-4,-4,NOBASE,NOBASE, 3,
    3,NOBASE,NOBASE,-4,NOBASE,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4 };
  if (d->ps >= d->psmax)
   if (d->psmax > 0L)
    { fseek(d->f,d->seqstart,SEEK_SET);
      d->ps = 0L; }
  NL:
  if ((ic = getc(d->f)) == EOF) goto FAIL;
  ic = map[ic];
  if (ic >= 0)
   { d->ps++;
     return(ic); }
  if (ic == -2)
   { d->nextseq = ftell(d->f) - 1L;
     return(TERM); }
  goto NL;
  FAIL:
  d->nextseq = -1L;
  if (d->psmax > 0L) 
   { d->ps = d->psmax;
     return(NOBASE); }
  else return(TERM); }


int seq_init(data_set *d)
{ long ngc;
  int ic;
  if ((d->seqstart = get_seqname(d)) < 0L) return(0);
  d->ps = 0L;
  d->psmax = -1L;
  ngc = 0L;
  while ((ic = move_forward(d)) >= 0) if (ic > 0) if (ic < 3) ngc++;
  if ((d->psmax = d->ps) <= 0L) return(0);
  d->gc = (double)ngc/(double)d->psmax;
  fseek(d->f,d->seqstart,SEEK_SET);
  d->ps = 0L;
  return(1); }



char cbase(int c)
{ static char base[6] = "acgt.^";
  if (c < 0) return('#');
  if (c > 5) return((char)c);
  return(base[c]); }


char cpbase(int c)
{ static char base[6] = "ACGT.^";
  if (c < 0) return('#');
  if (c > 5) return((char)c);
  return(base[c]); }




char *aa(int *anticodon, csw *sw)
{ int p1,p2,p3;
  static char aan[65][5] =
   { "Phe","Val","Leu","Ile",
     "Cys","Gly","Arg","Ser",
     "Ser","Ala","Pro","Thr",
     "Tyr","Asp","His","Asn",
     "Leu","Val","Leu","Met",
     "Trp","Gly","Arg","Arg",
     "Ser","Ala","Pro","Thr",
     "Stop","Glu","Gln","Lys",
     "Phe","Val","Leu","Ile",
     "Cys","Gly","Arg","Ser",
     "Ser","Ala","Pro","Thr",
     "Tyr","Asp","His","Asn",
     "Leu","Val","Leu","Ile",
     "Sec","Gly","Arg","Arg",
     "Ser","Ala","Pro","Thr",
     "Stop","Glu","Gln","Lys","" };
  static char mt_aan[65][5] =
   { "Phe","Val","Leu","Ile",
     "Cys","Gly","Arg","Ser",
     "Ser","Ala","Pro","Thr",
     "Tyr","Asp","His","Asn",
     "Leu","Val","Leu","Met",
     "Trp","Gly","Arg","Arg",
     "Ser","Ala","Pro","Thr",
     "Stop","Glu","Gln","Lys",
     "Phe","Val","Leu","Ile",
     "Cys","Gly","Arg","Ser",
     "Ser","Ala","Pro","Thr",
     "Tyr","Asp","His","Asn",
     "Leu","Val","Leu","Ile",
     "Trp","Gly","Arg","Arg",
     "Ser","Ala","Pro","Thr",
     "Stop","Glu","Gln","Lys","" };
  if ((p1 = *anticodon) >= NOBASE) return(aan[64]);
  if ((p2 = anticodon[1]) >= NOBASE) return(aan[64]);
  if ((p3 = anticodon[2]) >= NOBASE) return(aan[64]);
  if (sw->mt)
   return(mt_aan[(p1<<4) + (p2<<2) + p3]);
  else
   return(aan[(p1<<4) + (p2<<2) + p3]); }


char *translate(int *codon)
{ int p1,p2,p3;
  static char aan[65][5] =
   { "Phe","Val","Leu","Ile",
     "Cys","Gly","Arg","Ser",
     "Ser","Ala","Pro","Thr",
     "Tyr","Asp","His","Asn",
     "Leu","Val","Leu","Met",
     "Trp","Gly","Arg","Arg",
     "Ser","Ala","Pro","Thr",
     "Stop","Glu","Gln","Lys",
     "Phe","Val","Leu","Ile",
     "Cys","Gly","Arg","Ser",
     "Ser","Ala","Pro","Thr",
     "Tyr","Asp","His","Asn",
     "Leu","Val","Leu","Ile",
     "Stop","Gly","Arg","Arg",
     "Ser","Ala","Pro","Thr",
     "Stop","Glu","Gln","Lys",
     "" };
  if ((p1 = *codon) >= NOBASE) return(aan[64]);
  if ((p2 = codon[1]) >= NOBASE) return(aan[64]);
  if ((p3 = codon[2]) >= NOBASE) return(aan[64]);
  return(aan[((3 - p3)<<4) + ((3 - p2)<<2) + (3 - p1)]); }

char ltranslate(int *codon)
{ int p1,p2,p3;
  static char aan[65] =
   { "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF-" };
  if ((p1 = *codon) >= NOBASE) return(aan[64]);
  if ((p2 = codon[1]) >= NOBASE) return(aan[64]);
  if ((p3 = codon[2]) >= NOBASE) return(aan[64]);
  return(aan[(p1<<4) + (p2<<2) + p3]); }

char ptranslate(int *codon)
{ int p1,p2,p3;
  static char polar[105] =
   { "PPPPPPPPPPPPNNNNPPPPNNNNPPPPNNNNPPPPNNNNNNNNNNNN*N*NPPPP*NNNNNNN-" };
  if ((p1 = *codon) >= NOBASE) return(polar[64]);
  if ((p2 = codon[1]) >= NOBASE) return(polar[64]);
  if ((p3 = codon[2]) >= NOBASE) return(polar[64]);
  return(polar[(p1<<4) + (p2<<2) + p3]); }



double gc_content(gene *t)
{ int *s,*se;
  double ngc;
  static double score[6] = { 0.0,1.0,1.0,0.0,0.0,0.0 };
  ngc = 0.0;
  if ((t->nintron > 0) && (t->asst == 0))
   { s = t->eseq + ASTEM1_EXT;
     se = s + t->intron;
     while (s < se) ngc += score[*s++];
     s = se + t->nintron;
     se = t->eseq + ASTEM1_EXT + t->nbase + t->nintron;
     while (s < se) ngc += score[*s++]; }
  else
   { s = (t->status == 3)?t->eseq:t->seq;
     se = s + t->nbase;
     while (s < se) ngc += score[*s++]; }
  return(ngc/(double)t->nbase); }


void init_tmrna(FILE *f, csw *sw)
{ int i,c,*s;
  s = sw->tmrna_struct;
  while ((c = *s++) != TERM) itmparam(cbase(c),f); }


int *make_var(int *seq, char matrix[][MATY],
               int *x, int *y, int orient, int var)
{ int i,px,py,l,stem;
  static int ux[4] = { 1,0,-1,0 };
  static int uy[4] = { 0,1,0,-1 };
  static int vx[4] = { 0,-1,0,1 };
  static int vy[4] = { 1,0,-1,0 };
  static int loopu[8][8] =
  { { 2 }, { 1,1 }, { 1,1,0 }, { 1,0,1,0 }, { 1,1,0,0,0 },
    { 1,1,1,-1,-1,1 }, { 1,1,1,0,0,-1,0 }, { 1,1,1,1,0,-1,-1,0 } };
  static int loopv[8][8] =
  { { 3 }, { 1,2 },  { 0,1,2 }, { 0,1,1,1 }, { -1,1,1,1,1 },
    { -1,0,1,1,1,1 }, { -1,-1,1,1,1,1,1 }, { -1,-1,0,1,1,1,1,1 } };
  px = *x;
  py = *y;
  if (var < 0) var = 0;
  if (var > 30) var = 30;
  if (var % 2)
   stem = (var - 7)/2;
  else
   stem = (var - 6)/2;
  if (stem < 0) stem = 0;
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
  return(seq); }


int *make_tv(int *seq, char matrix[][MATY],
             int *x, int *y, int orient, int tv)
{ int i,px,py,l,stem;
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
  static char base1[7] = "acgtgt";
  static char base2[7] = "tgcatg";
  static int score[7] = { 2,2,2,2,1,1 };
  s = 0;
  for (i = 0; i < 6; i++)
   if (b1 == base1[i])
    if (b2 == base2[i])
     { s = score[i];
       break; }
  return(s); }


int *make_clover(int *seq, int b, int e, int stemlength,
                  char matrix[][MATY], int *x, int *y, int orient)
{ int i,px,py,l;
  int *s,*se;
  static int ux[5] = { 1,0,-1,0,0 };
  static int uy[5] = { 0,1,0,-1,1 };
  static int vx[5] = { 0,-1,0,1,1 };
  static int vy[5] = { 1,0,-1,0,0 };
  static int loopu[17][17] =
  { { -1 }, { 0,-1 }, { 0,0,-1 }, { 0,1,-1,-1 }, { 0,1,0,-1,-1 },
    { 0,1,0,0,-1,-1 }, { 0,1,1,0,-1,-1,-1 }, { 0,1,1,0,0,-1,-1,-1 },
    { 0,1,1,1,0,-1,-1,-1,-1 }, { 0,1,1,1,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,0,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 } };
  static int loopv[17][17] =
  { { 2 }, { 1,1 },  { 0,1,1 }, { -1,2,2,-1 }, { -1,1,1,2,-1 },
    { -1,1,1,1,1,-1 }, { -1,0,1,1,1,1,-1 }, { -1,0,1,1,1,1,0,-1 },
    { -1,0,1,1,1,1,0,0,-1 }, { -1,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 } };
  static int dloopu[17][17] =
  { { -1 }, { 0,-1 }, { 0,0,-1 }, { 0,1,-1,-1 }, { 0,1,0,-1,-1 },
    { 0,1,0,0,-1,-1 }, { 0,1,1,0,-1,-1,-1 }, { 0,1,1,0,0,-1,-1,-1 },
    { 0,1,1,0,0,0,-1,-1,-1 }, { 0,1,1,1,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,0,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 } };
  static int dloopv[17][17] =
  { { 2 }, { 1,1 },  { 0,1,1 }, { -1,2,2,-1 }, { -1,1,1,2,-1 },
    { -1,1,1,1,1,-1 }, { -1,0,1,1,1,1,-1 }, { -1,0,1,1,1,1,0,-1 },
    { -1,0,1,1,1,1,1,-1,-1 }, { -1,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 } };
  static char bond1[4] = " +!";
  static char bond2[4] = " +-";
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
  if (l < 17)
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
   { px += ux[orient]*loopu[0][0] + vx[orient]*loopv[0][0];
     py += uy[orient]*loopu[0][0] + vy[orient]*loopv[0][0]; }
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
       case  UPRIGHT:
       case  UP:
       case  DOWN:  matrix[px - vx[orient]][py - vy[orient]] = bond2[i];
                    break; }
    px -= ux[orient];
    py -= uy[orient]; }
  *x = px;
  *y = py;
  return(seq+e); }





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



void build_trna(gene *t, char matrix[][MATY], int x, int y)
{ int i,j,e,c,*seq; /* *anticodon */
  int rseq[150];
  static char bond2[4] = " +-";
  remove_inserts(t->seq,rseq);
  seq = rseq;
  i = 0;
  while (i < t->astem1)
  { matrix[x][y] = cbase(*seq++);
    y--;
    i++; }
  if (t->spacer > 0)
   { x--;
     matrix[x][y] = cbase(*seq++);
     x--;
     if (t->spacer == 3) matrix[x][y] = cbase(*seq++);
     y--;
     if (t->spacer >= 2) matrix[x][y] = cbase(*seq++);
     x--;
     y--; }
  if (t->dstem > 0)
   { e = 2*t->dstem + t->dloop;
     seq = make_clover(seq,0,e,t->dstem,matrix,&x,&y,LEFT);
     y--;
     matrix[x][y] = cbase(*seq++);
     x++;
     y--; }
  else
   seq = make_dv(seq,matrix,t->dloop,RIGHT,&x,&y);
  e = 2*t->cstem + t->cloop;
  /* anticodon = seq + t->cstem + 2; */
  seq = make_clover(seq,0,e,t->cstem,matrix,&x,&y,DOWN);
  if (t->tstem > 0)
   { seq = make_var(seq,matrix,&x,&y,RIGHT,t->var);
     e = 2*t->tstem + t->tloop;
     seq = make_clover(seq,0,e,t->tstem,matrix,&x,&y,RIGHT);
     y++; }
  else
   seq = make_tv(seq,matrix,&x,&y,RIGHT,t->tloop);
  e = t->astem2;
  i = 0;
  while (i < e)
  { if ((c = *seq++) < 0) break;
    matrix[x][y] = cbase(c);
    j = base_match(matrix[x][y],matrix[x - 2][y]);
    matrix[x - 1][y] = bond2[j];
    y++;
    i++; }
  i = 0;
  while (i < 2)
  { if ((c = *seq++) < 0) break;
    matrix[x][y] = cbase(c);
    x++;
    y++;
    i++; }
  e = ASTEM2_EXTD-2;
  i = 0;
  while (i < e)
  { if ((c = *seq++) < 0) break;
    matrix[x][y] = cbase(c);
    x++;
    i++; }
  /* return(anticodon); */ }




void build_tmrna(gene *t, char matrix[][MATY], int x, int y)
{ int i,j,e,c,tarm,*seq;
  int rseq[2*MAXTMRNALEN+1];
  static char bond2[4] = " +-";
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
  seq = make_var(seq,matrix,&x,&y,RIGHT,t->var);
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
  i = 0;
  while (i < 2)
  { if ((c = *seq++) == TERM) break;
    matrix[x][y] = cbase(c);
    x++;
    y++;
    i++; }
  e = ASTEM2_EXTD-2;
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
{ int i,j;
  i = ylines;
  while (--i >= 0)
   { for (j = 0; j < MATX; j++) fputc(matrix[j][i],f);
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
  char *s,*st,*sb,*sd,*sdb;
  static struct { char name[50]; char tag[50]; } tagdatabase[NTAG] =
   { { "Cyanidioschyzon merolae Chloroplast","ANQILPFSIPVKHLAV" },
     { "Mesostigma viride chloroplast","ANNILPFNRKTAVAV" },
     { "Nephroselmis olivacea chloroplast","TTYHSCLEGHLS" },
     { "Pirellula sp.","AEENFALAA" },
     { "Desulfotalea psychrophila","ADDYNYAVAA" },
     { "Desulfuromonas acetoxidans","ADTDVSYALAA" },
     { "Exiguobacterium sp.","GKTNTQLAAA" },
     { "Mycoplasma gallisepticum","DKTSKELADENFVLNQLASNNYALNF" },
     { "Aquifex aeolicus","APEAELALAA" },
     { "Thermotoga maritima ","ANEPVAVAA" },
     { "Chloroflexus aurantiacus","ANTNTRAQARLALAA" },
     { "Thermus thermophilus","ANTNYALAA" },
     { "Deinococcus radiodurans","GNQNYALAA" },
     { "Cytophaga hutchinsonii","GEESYAMAA" },
     { "Bacteroides fragilis","GETNYALAA" },
     { "Tannerella forsythensis","GENNYALAA" },
     { "Porphyromonas gingivalis","GENNYALAA" },
     { "Prevotella intermedia","GENNYALAA" },
     { "Chlorobium tepidum","ADDYSYAMAA" },
     { "Gemmata obscuriglobus","AEPQYSLAA" },
     { "Chlammydophila pneumoniae","AEPKAECEIISLFDSVEERLAA" },
     { "Chlammydophila caviae","AEPKAECEIISFSDLTEERLAA" },
     { "Chlammydophila abortus","AEPKAKCEIISFSELSEQRLAA" },
     { "Chlammydia trachomatis","AEPKAECEIISFADLEDLRVAA" },
     { "Chlammydia muridarum","AEPKAECEIISFADLNDLRVAA" },
     { "Nostoc PCC7120","ANNIVKFARKDALVAA" },
     { "Nostoc punctiforme","ANNIVNFARKDALVAA" },
     { "Fremyella diplosiphon","ANNIVKFARKEALVAA" },
     { "Plectonema boryanum","ANNIVPFARKTAPVAA" },
     { "Trichodesmium erythraeum","ANNIVPFARKQVAALA" },
     { "Oscillatoria 6304","ANNIVPFARKAAPVAA" },
     { "Chroococcidiopsis PCC6712","ANNIVKFERQAVFA" },
     { "Synechocystis 6803","ANNIVSFKRVAIAA" },
     { "Thermosynechococcus elongatus","ANNIVPFARKAAAVA" },
     { "Synechococcus PCC6301","ANNIVPFARKAAPVAA" },
     { "Synechococcus WH8102","ANNIVRFSRHAAPVAA" },
     { "Synechococcus PCC6307","ANNIVRFSRQAAPVAA" },
     { "Synechococcus PCC7002","ANNIVPFARKAAAVA" },
     { "Synechococcus PCC7009","ANNIVRFSRQAAPVAA" },
     { "Synechococcus PCC6904","ANNIVRFSRQAAPVAA" },
     { "Prochlorococcus marinus 1","ANKIVSFSRQTAPVAA" },
     { "Prochlorococcus marinus 2","ANNIVRFSRQPALVAA" },
     { "Prochlorococcus marinus 3","ANKIVSFSRQTAPVAA" },
     { "Cyanophora paradoxa chloroplast","ATNIVRFNRKAAFAV" },
     { "Thalassiosira weissflogii chloroplast","ANNIIPFIFKAVKTKKEAMALNFAV" },
     { "Odontella sinensis chloroplast","ANNLISSVFKSLSTKQNSLNLSFAV" },
     { "Bolidomonas pacifica chloroplast","ANNILAFNRKSLSFA" },
     { "Pavlova lutheri chloroplast","ANNILSFNRVAVA" },
     { "Porphyra purpurea chloroplast","AENNIIAFSRKLAVA" },
     { "Guillardia theta chloroplast","ASNIVSFSSKRLVSFA" },
     { "Fibrobacter succinogenes","ADENYALAA" },
     { "Treponema pallidum","ANSDSFDYALAA" },
     { "Treponema denticola","AENNDSFDYALAA" },
     { "Leptospira interrogans","ANNELALAA" },
     { "Borrelia burgdorferi","AKNNNFTSSNLVMAA" },
     { "Caulobacter crescentus","ANDNFAEEFAVAA" },
     { "Rhodobacter sphaeroides","ANDNRAPVALAA" },
     { "Silicibacter pomeroyi","ANDNRAPVALAA" },
     { "Rhodopseudomonas palustris","ANDNYAPVAQAA" },
     { "Bradyrhizobium japonicum","ANDNFAPVAQAA" },
     { "Agrobacterium tumefaciens 1","ANDNNAKEYALAA" },
     { "Agrobacterium tumefaciens 2","ANDNNAKECALAA" },
     { "Rhizobium leguminosarum","ANDNYAEARLAA" },
     { "Sinorhizobium meliloti","ANDNYAEARLAA" },
     { "Mesorhizobium loti","ANDNYAEARLAA" },
     { "Brucella melitensis","ANDNNAQGYALAA" },
     { "Brucella abortus","ANDNNAQGYALAA" },
     { "Methylobacterium extorquens","ANDNFAPVAVAA" },
     { "Magnetospirillum magnetotacticum","ANDNVELAAAA" },
     { "Rhodospirillum rubrum","ANDNVELAAAA" },
     { "Novosphingobium aromaticivorans","ANDNEALALAA" },
     { "Ehrlichia chaffeensis","ANDNFVFANDNNSSANLVAA" },
     { "Anaplasma phagocytophilum","ANDDFVAANDNVETAFVAAA" },
     { "Wolbachi.sp","ANDNFAAEDNVDAIAA" },
     { "Rickettsia conorii","ANDNNRSVGHLALAA" },
     { "Rickettsia typhi","ANDNKRYVGVAALAAA" },
     { "Rickettsia prowazekii","ANDNRYVGVPALAAA" },
     { "Neisseria gonorrhoeae","ANDETYALAA" },
     { "Neisseria meningitidis","ANDETYALAA" },
     { "Chromobacterium violaceum","ANDETYALAA" },
     { "Uncultured U02","ANDEQFALAA" },
     { "Nitrosomonas europaea","ANDENYALAA" },
     { "Nitrosomonas cryotolerans","ANDENYALAA" },
     { "Methylobacillus glycogenes","ANDETYALAA" },
     { "Uncultured U04","ANDETYALAA" },
     { "Ralstonia pickettii","ANDERYALAA" },
     { "Ralstonia solanacearum","ANDNRYQLAA" },
     { "Ralstonia eutropha","ANDERYALAA" },
     { "Ralstonia metallidurans","ANDERYALAA" },
     { "Alcaligenes faecalis","ANDERFALAA" },
     { "Comamonas testosteroni","ANDERFALAA" },
     { "Variovorax paradoxus","ANDERFALAA" },
     { "Hydrogenophaga palleronii","ANDERFALAA" },
     { "Burkholderia pseudomallei","ANDDTFALAA" },
     { "Burkholderia mallei","ANDDTFALAA" },
     { "Burkholderia fungorum","ANDDTFALAA" },
     { "Burkholderia cepacia","ANDDTFALAA" },
     { "Bordetella pertussis","ANDERFALAA" },
     { "Bordetella parapertussis","ANDERFALAA" },
     { "Bordetella bronchiseptica","ANDERFALAA" },
     { "Xylella fastidiosa 1","ANEDNFAVAA" },
     { "Xylella fastidiosa 2","ANEDNFALAA" },
     { "Xylella fastidiosa 3","ANEDNFAIAA" },
     { "Xanthomonas campestris","ANDDNYGSDFAIAA" },
     { "Xanthomonas axonopodis","ANDDNYGSDFAIAA" },
     { "Legionella pneumophila","ANDENFAGGEAIAA" },
     { "Coxiella burnetii","ANDSNYLQEAYA" },
     { "Methylococcus capsulatus","ANDDVYALAA" },
     { "Uncultured U01a","ANDSNYALAA" },
     { "Dichelobacter nodosus","ANDDNYALAA" },
     { "Francisella tularensis","GNKKANRVAANDSNFAAVAKAA" },
     { "Acidithiobacillus ferrooxidans","ANDSNYALAA" },
     { "Acinetobacter ADP1","ANDETYALAA" },
     { "Psychrobacter 2734","ANDENYALAA" },
     { "Azotobacter vinelandii","ANDDNYALAA" },
     { "Pseudomonas aeruginosa","ANDDNYALAA" },
     { "Pseudomonas syringae 1","ANDENYGAQLAA" },
     { "Pseudomonas syringae 2","ANDETYGEYALAA" },
     { "Pseudomonas fluorescens 1","ANDDQYGAALAA" },
     { "Pseudomonas fluorescens 2","ANDENYGQEFALAA" },
     { "Pseudomonas putida 1","ANDENYGAEYKLAA" },
     { "Marinobacter hydrocarbonoclasticus","ANDENYALAA" },
     { "Pseudoalteromonas haloplanktis","ANDDNYSLAA" },
     { "Uncultured WW11","ANDDNYALAA" },
     { "Shewanella oneidensis","ANDDNYALAA" },
     { "Photorhabdus asymbiotica?","ANDNEYALVA" },
     { "Microbulbifer degradans","ANDDNYGAQLAA" },
     { "Colwellia sp","ANDDTFALAA" },
     { "Photobacterium phosphoreum","ANDENYALAA" },
     { "Vibrio cholerae","ANDENYALAA" },
     { "Vibrio vulnificus","ANDENYALAA" },
     { "Aeromonas salmonicida","ANDENYALAA" },
     { "Uncultured VLW3","ANDENYALAA" },
     { "Uncultured VLS13","ANDENYALAA" },
     { "Uncultured WW9","ANDENYALAA" },
     { "Uncultured WW10","ANDENYALAV" },
     { "Uncultured VLW5","ANDENYALAA" },
     { "Uncultured RCA4","ANDETYALAA" },
     { "Uncultured LEM1","ANDETYALAA" },
     { "Uncultured LEM2","ANDETHALAA" },
     { "Wigglesworthia brevipalpis","AKHKYNEPALLAA" },
     { "Buchnera aphidicola 1","ANNKQNYALAA" },
     { "Buchnera aphidicola 2","ANNKQNYALAA" },
     { "Buchnera aphidicola 3","AKQNQYALAA" },
     { "Shigella dysenteriae","ANDENYALAA" },
     { "Shigella flexneri","ANDENYALAA" },
     { "Escherichia coli","ANDENYALAA" },
     { "Providencia rettgeri","ANDENYALAA" },
     { "Serratia marcescens","ANDENYALAA" },
     { "Klebsiella pneumoniae","ANDENYALAA" },
     { "Pectobacterium carotovora","ANDENYALAA" },
     { "Erwinia chrysanthemi","ANDENFAPAALAA" },
     { "Salmonella bongori","ANDENYALAA" },
     { "Salmonella typhimurium","ANDETYALAA" },
     { "Salmonella typhi","ANDETYALAA" },
     { "Salmonella paratyphi","ANDENYALAA" },
     { "Salmonella enterica 1","ANDETYALAA" },
     { "Salmonella enterica 2","ANDENYALAA" },
     { "Uncultured RCA1","ANDENYALAA" },
     { "Uncultured VLS1","ANDENYALAA" },
     { "Uncultured WW1","ANDENYALAA" },
     { "Uncultured RCA2","SNDENYALAA" },
     { "Uncultured WW2","ANDENYALAA" },
     { "Uncultured QL1","ANVENYALAA" },
     { "Uncultured WW4","ANDGNYALAA" },
     { "Uncultured VLS5","ANDETYALAA" },
     { "Uncultured FS1","ANDETYALAA" },
     { "Uncultured VLS6","ANDENYALAA" },
     { "Uncultured FS2","ANDENYALAA" },
     { "Uncultured WW5","ANDENYALAA" },
     { "Uncultured VLW1","ANDENYALAA" },
     { "Uncultured VLS7","ANDENYALAA" },
     { "Uncultured VLS9","ANDENYALAA" },
     { "Uncultured VLW2","ANDENYALAA" },
     { "Uncultured WW7","ANDENCALAA" },
     { "Uncultured WW8","ANDENYALAA" },
     { "Yersinia pestis","ANDENYALAA" },
     { "Yersinia enterocolitica","ANDSQYESAALAA" },
     { "Mannheimia haemolytica","ANDEQYALAA" },
     { "Haemophilus ducreyi","ANDEQYALAA" },
     { "Haemophilus influenzae","ANDEQYALAA" },
     { "Haemophilus somnus","ANDEQYALAA" },
     { "Pasteurella multocida","ANDEQYALAA" },
     { "Actinobacillus actinomycetemcomitans","ANDEQYALAA" },
     { "Actinobacillus pleuropneumoniae","ANDEQYALAA" },
     { "Lawsonia intracellularis","ANNNYDYALAA" },
     { "Desulfovibrio desulfuricans","ANNDYDYAYAA" },
     { "Desulfovibrio vulgaris","ANNYDYALAA" },
     { "Geobacter sulfurreducens","ADNYDYAVAA" },
     { "Geobacter metallireducens","ADNYDYAVAA" },
     { "Helicobacter pylori 1","VNNTDYAPAYAKAA" },
     { "Helicobacter pylori 2","VNNADYAPAYAKAA" },
     { "Helicobacter pylori 3","VNNTDYAPAYAKAA" },
     { "Campylobacter jejuni","ANNVKFAPAYAKAA" },
     { "Fusobacterium nucleatum 1","GNKDYALAA" },
     { "Fusobacterium nucleatum 2","GNKEYALAA" },
     { "Dehalococcoides ethenogenes","GERELVLAG" },
     { "Mycobacterium leprae","ADSYQRDYALAA" },
     { "Mycobacterium avium","ADSHQRDYALAA" },
     { "Mycobacterium bovis","ADSHQRDYALAA" },
     { "Mycobacterium tuberculosis","ADSHQRDYALAA" },
     { "Mycobacterium marinum","ADSHQRDYALAA" },
     { "Mycobacterium smegmatis","ADSNQRDYALAA" },
     { "Corynebacterium diphtheriae","AENTQRDYALAA" },
     { "Corynebacterium glutamicum","AEKSQRDYALAA" },
     { "Thermobifida fusca","ANSKRTEFALAA" },
     { "Streptomyces coelicolor","ANTKRDSSQQAFALAA" },
     { "Tropheryma whipplei","ANLKRTDLSLAA" },
     { "Clavibacter michiganensis","ANNKQSSFVLAA" },
     { "Bifidobacterium longum","AKSNRTEFALAA" },
     { "Bifidobacterium longum","AKSNRTEFALAA" },
     { "Bacillus anthracis","GKQNNLSLAA" },
     { "Bacillus cereus","GKQNNLSLAA" },
     { "Bacillus megaterium","GKSNNNFALAA" },
     { "Bacillus halodurans","GKENNNFALAA" },
     { "Bacillus subtilis","GKTNSFNQNVALAA" },
     { "Bacillus stearothermophilus","GKQNYALAA" },
     { "Staphylococcus aureus","GKSNNNFAVAA" },
     { "Staphylococcus saprophyticus","GKENNNFAVAA" },
     { "Staphylococcus xylosus","GKENNNFAVAA" },
     { "Staphylococcus epidermidis","DKSNNNFAVAA" },
     { "Oceanobacillus iheyensis","GKETNQPVLAAA" },
     { "Listeria monocytogenes","GKEKQNLAFAA" },
     { "Listeria innocua","GKEKQNLAFAA" },
     { "Listeria welshimeri","GKEKQNLAFAA" },
     { "Listeria seeligeri","GKEKQNLAFAA" },
     { "Listeria grayi 1","GKEKQNLAFAA" },
     { "Listeria grayi 2","GKQNNNLAFAA" },
     { "Listeria ivanovii","GKEKQNLAFAA" },
     { "Lactobacillus gasseri","ANNENSYAVAA" },
     { "Lactobacillus sakei","ANNNNSYAVAA" },
     { "Lactobacillus helveticus","ANNKNSYALAA" },
     { "Lactobacillus gallinarum","ANNKNSYALAA" },
     { "Lactobacillus plantarum","AKNNNNSYALAA" },
     { "Pediococcus pentosaceus","AKNNNNSYALAA" },
     { "Leuconostoc mesenteroides","AKNENSFAIAA" },
     { "Leuconostoc lactis","AKNENSFAIAA" },
     { "Leuconostoc pseudomesenteroides","AKNENSYAIAA" },
     { "Enterococcus durans","AKNENNSYALAA" },
     { "Oenococcus oeni","AKNNEPSYALAA" },
     { "Enterococcus faecium","AKNENNSYALAA" },
     { "Enterococcus faecalis","AKNENNSFALAA" },
     { "Streptococcus equi","AKNNTTYALAA" },
     { "Streptococcus suis","AKNTNTYALAA" },
     { "Streptococcus uberis","AKNTNSYALAA" },
     { "Streptococcus pyogenes","AKNTNSYALAA" },
     { "Streptococcus agalactiae","AKNTNSYALAA" },
     { "Streptococcus mutans","AKNTNSYAVAA" },
     { "Streptococcus gordonii","AKNNTSYALAA" },
     { "Streptococcus pneumoniae","AKNNTSYALAA" },
     { "Streptococcus mitis","AKNNTSYALAA" },
     { "Streptococcus thermophilus","AKNTNSYAVAA" },
     { "Lactococcus raffinolactis","AKNTQTYAVAA" },
     { "Lactococcus plantarum","AKNTQTYALAA" },
     { "Lactococcus garvieae","AKNNTSYALAA" },
     { "Lactococcus lactis","AKNNTQTYAMAA" },
     { "Mycoplasma capricolum","ANKNEETFEMPAFMMNNASAGANFMFA" },
     { "Mesoplasma florum","ANKNEENTNEVPTFMLNAGQANYAFA" },
     { "Spiroplasma kunkelii","ASKKQKEDKIEMPAFMMNNQLAVSMLAA" },
     { "Ureaplasma urealyticum","AENKKSSEVELNPAFMASATNANYAFAY" },
     { "Mycoplasma pulmonis","GTKKQENDYQDLMISQNLNQNLAFASV" },
     { "Mycoplasma penetrans","AKNNKNEAVEVELNDFEINALSQNANLALYA" },
     { "Mycoplasma genitalium 1","DKENNEVLVEPNLIINQQASVNFAFA" },
     { "Mycoplasma genitalium 2","DKENNEVLVDPNLIINQQASVNFAFA" },
     { "Mycoplasma pneumoniae","DKNNDEVLVDPMLIANQQASINYAFA" },
     { "Thermoanaerobacter tengcongensis","ADRELAYAA" },
     { "Heliobacillus mobilis","AEDNYALAA" },
     { "Desulfitobacterium hafniense","ANDDNYALAA" },
     { "Carboxydothermus hydrogenoformans","ANENYALAA" },
     { "Ruminococcus albus","GHGYFAKAS" },
     { "Clostridium acetobutylicum","DNENNLALAA" },
     { "Clostridium perfringens","AEDNFALAA" },
     { "Clostridium thermocellum","ANEDNYALAAA" },
     { "Clostridium botulinum","ANDNFALAA" },
     { "Clostridium tetani","ADDNFVLAA" },
     { "Clostridium difficile","ADDNFAIAA" },
     { "Hyphomonas neptunium","ANDNFAEGELLAA" },
     { "Vibrio fischeri","ANDENYALAA" },
     { "Corynebacterium efficiens","AEKTQRDYALAA" },
     { "Streptomyces avermitilus","ANTKSDSQSFALAA" },
     { "Brevibacterium linens","AKSNNRTDFALAA" },
     { "Lactobacillus delbrueckii 1","AKNENNSYALAA" },
     { "Lactobacillus delbrueckii 2","ANENSYAVAA" },
     { "Lactobacillus casei","AKNENSYALAA" },
     { "Lactobacillus brevis","AKNNNNSYALAA" },
     { "Streptomyces thermophilus","AKNTNSYAVAA" } };
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
  copy("many bacteria",thit[0]);
  return(1);  }
            

void disp_peptide_tag(FILE *f, gene *t, csw *sw)
{ int i,j,lx,nm,nmh,c1,c2,c3,*s,*se;
  char tag[50],thit[21][50];
  fprintf(f,"Tag peptide (at %d)\nTag sequence: ",t->tps+1);
  se = t->eseq + t->tps;
  lx = (t->tpe - t->tps + 1);
  if (ltranslate(se+lx) == '*')
   { lx += 3;
     if (ltranslate(se+lx) == '*') lx += 3; }
  lx /= 3;
  s = se;
  for (i = 0; i < lx; i++)
  { if (i > 0) fputc('-',f);
    if ((c1 = *s++) >= NOBASE) continue; 
    if ((c2 = *s++) >= NOBASE) continue; 
    if ((c3 = *s++) >= NOBASE) continue; 
    fputc(cbase(c1),f);
    fputc(cbase(c2),f);
    fputc(cbase(c3),f); }
  s = se;
  fprintf(f,"\nTag peptide:  ");
  for (i = 0; i < lx; i++)
  { fprintf(f,"%s",translate(s));
    s += 3;
    if (i < (lx-1)) fputc('-',f); }
  s = se;
  fprintf(f,"\nTag peptide:  ");
  for (i = 0; i < lx; i++)
  { tag[i] = ltranslate(s);
    fprintf(f,"%c",tag[i]);
    s += 3; }
  tag[lx] = '\0';
  if (sw->energydisp)
   { s = se;
     fprintf(f,"\nTag Polarity: ");
     for (i = 0; i < lx; i++)
     { fprintf(f,"%c",ptranslate(s));
       s += 3; }}
  fputc('\n',f);
  nmh = identify_tag(tag,lx,thit,21);
  if (nmh > 0)
   { if (nmh > 1)
      { fprintf(f,"Tag match to tmRNA's from:\n");
        i = 0;
        for (nm = 0; nm < nmh; nm++)
         { if (++i > 3) 
            { fputc('\n',f);
              i = 1; }
            if (i > 1) fprintf(f,", ");
           fprintf(f,thit[nm]); }
        fputc('\n',f); }
     else
       fprintf(f,"Tag match to tmRNA from %s\n",thit[0]); }
  else fprintf(f,"Tag not identified\n");
  fputc('\n',f);  }


void sense_switch(int *seq1, int *seq2, int lseq)
{ int i;
  static int cmap[6] = { 3,2,1,0,4,5 };
  seq2 += lseq;
  i = -1;
  while (++i < lseq) *--seq2 = cmap[*seq1++]; }
	

void disp_intron(FILE *f, gene *t)
{ int i,k,p,c,no,noh,nr,nrh,*s,*sb,*se,*sf,*so;
  int cseq[MAXINTRONLEN];
  astem ohit[50];
  trna_loop rhit[50];
  if (t->nintron <= 0) return;
  fprintf(f,"Intron from %s\n",t->name);
  fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n");
  sb = t->eseq + ASTEM1_EXT + t->intron;
  s = sb;
  se = sb + t->nintron;
  i = 0;
  while (s < se)
   { if ((c = *s++) < 0) break;
     fputc(cbase(c),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f);
  fputc('\n',f);
  fprintf(f,"Intron Length: %d\n",t->nintron);
  fprintf(f,"Intron Insertion Position(%d): ",t->intron+1);
  s = sb - 5;
  for (i = 0; i < 5; i++) fputc(cbase(*s++),f);
  fprintf(f,"-Intron-");
  s = se;
  for (i = 0; i < 5; i++) fputc(cbase(*s++),f);
  fputc('\n',f);
  fputc('\n',f); }


void disp_seq(FILE *f, gene *t, csw *sw)
{ int i,*s,*se;
  if (!sw->batch)
   fprintf(f,"\n1   .   10    .   20    .   30    .   40    .   50\n");
  if (t->nintron > 0)
   { s = t->eseq;
     se = s + ASTEM1_EXT + t->nbase + t->nintron; }
  else
   { s = t->seq;
     se = s + ASTEM1_EXT + t->nbase; }
  i = 0;
  while (s < se)
   { fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f);
  if (!sw->batch) fputc('\n',f); }



void disp_tmrna_seq(FILE *f, gene *t, csw *sw)
{ int i,j,k,lx,c,*s,*sb,*se,*cseq;
  if (t->nintron <= 0) return;
  if (*(t->name) == '\0') fprintf(f,"tmRNA Sequence\n\n");
  else fprintf(f,"tmRNA Sequence in %s\n\n",t->name);
  fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n");
  sb = t->eseq + ASTEM1_EXT;
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
  while (ltranslate(se) == '*') se += 3;
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
  fputc('\n',f);
  fprintf(f,"Resume consensus sequence (at %d): ",t->tps - 6);
  s = t->eseq + t->tps - 7;
  for (i = 0; i < 18; i++) fputc(cbase(*s++),f);
  fputc('\n',f);
  fputc('\n',f);
  disp_peptide_tag(f,t,sw); }



void disp_tmrna_perm_seq(FILE *f, gene *t, csw *sw)
{ int i,j,k,lx,c,*s,*sb,*se,*cseq;
  if (t->nintron <= 0) return;
  if (*(t->name) == '\0') fprintf(f,"tmRNA Sequence\n\n");
  else fprintf(f,"tmRNA Sequence in %s\n\n",t->name);
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
  while (ltranslate(se) == '*') se += 3;
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
  fprintf(f,"\nResume consensus sequence (at %d): ",t->tps - 6);
  s = t->eseq + t->tps - 7;
  for (i = 0; i < 18; i++) fputc(cbase(*s++),f);
  fputc('\n',f);
  fputc('\n',f);
  disp_peptide_tag(f,t,sw); }

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
{ int *s;
  char *sb;
  switch (t->status)
   { case 3:  sprintf(si,"Stem-Loop");
              break;
     case 2:  sprintf(si,"Peptide Tag");
              break;
     case 1:  if (t->asst > 0)
               sprintf(si,"tmRNA (Permuted)");
              else
               sprintf(si,"tmRNA");
              break;
     case 0:  s = (proc?t->seq:t->ps) + t->anticodon;
              sb = si;
	      if (t->dstem == 0)
	       { sprintf(sb,"D-loop ");
		 sb += 7; }
	      if (t->tstem == 0)
	       { sprintf(sb,"TV-loop ");
		 sb += 8; }
              if (t->cloop == 8)
               sprintf(sb,"tRNA-%s(%c%c%c) or %s(%c%c%c)",
                          aa(s,sw),cbase(*s),cbase(s[1]),cbase(s[2]),
                          aa(s+1,sw),cbase(s[1]),cbase(s[2]),cbase(s[3]));
              else
               sprintf(sb,"tRNA-%s(%c%c%c)",
                          aa(s,sw),cbase(*s),cbase(s[1]),cbase(s[2]));
     default: break; }
  return(si); }



double vloop_stability(int *sv, int var)
{ int *sb,*se,nbp;
  unsigned int r,c,mx;
  static unsigned int A[6] = { 0,0,0,1,0,0 };
  static unsigned int C[6] = { 0,0,1,0,0,0 };
  static unsigned int G[6] = { 0,1,0,1,0,0 };
  static unsigned int T[6] = { 1,0,1,0,0,0 };
  unsigned int template[6];
  sb = sv + 3;
  se = sv + 7;
  template[0] = A[*sb];
  template[1] = C[*sb];
  template[2] = G[*sb];
  template[3] = T[*sb];
  template[4] = 0;
  template[5] = 0;
  while (++sb < se)
  { template[0] = (template[0] << 3) | A[*sb];
    template[1] = (template[1] << 3) | C[*sb];
    template[2] = (template[2] << 3) | G[*sb];
    template[3] = (template[3] << 3) | T[*sb]; }
  se +=  (((var - 4) >> 1) << 1) - 10;
  sb = se + 4;
  r = template[*sb];
  while (--sb > se) r = (r >> 3) + template[*sb];
  mx = r & 7;
  r = (r >> 3) + template[*sb--];
  c = r & 7;
  if (c > mx) mx = c;
  r = (r >> 3) + template[*sb];
  c = r & 7;
  if (c > mx) mx = c;
  return((double)(3*((int)mx - 4))); }



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



void disp_gene(gene *t, char m[][MATY], csw *sw)
{ int i,x,y,ng;
  double gc;
  char stat[80],*s;
  switch(t->status)
   { case 3: x = 13;
             y = 6;
             make_clover(t->eseq,0,2*t->tstem+t->tloop,
                         t->tstem,m,&x,&y,RIGHT);
             xcopy(m,11,3,"Stem-Loop",9);
             break; 
     case 2: x = 13;
             y = 6;
             make_clover(t->eseq,t->asst,t->asst+2*t->tstem+t->tloop,
                         t->tstem,m,&x,&y,RIGHT);
             xcopy(m,11,3,"Stem-Loop",9);
             break; 
     case 1: build_tmrna(t,m,13,27);
             xcopy(m,5,3,"tmRNA (tRNA domain)",19);
             break; 
     case 0: build_trna(t,m,13,27);
             name(t,stat,1,sw);
             xcopy(m,5,3,stat,length(stat));
             break; }
  location(stat,t,sw,"Sequence");
  xcopy(m,5,1,stat,length(stat));
  gc = gc_content(t);
  sprintf(stat,"%d bases, %%GC = %2.1lf",t->nbase,100.0*gc);
  xcopy(m,5,2,stat,length(stat));
  if (sw->energydisp)
   { sprintf(stat,"Score = %lg\n",t->energy);
     xcopy(m,5,0,stat,length(stat)); }}

void disp_batch_trna(FILE *f, gene *t, csw *sw)
{ int *s,anticodon;
  char pos[50];
  static char intron[3] = " I";
  position(pos,t,sw);
  anticodon = 1 + t->anticodon;
  if (t->nintron > 0)
   if (t->intron <= t->anticodon)
    anticodon += t->nintron;
  fprintf(f,"T%c%25s\t%d\t",intron[(t->nintron > 0)?1:0],pos,anticodon);
  s = t->seq + t->anticodon;
  fprintf(f,"%s(%c%c%c)",
             aa(s,sw),cbase(*s),cbase(s[1]),cbase(s[2]));
  if (t->nintron > 0)
   fprintf(f,"\t%d\t%d",t->intron+1,t->nintron);
  if (sw->energydisp)
   fprintf(f,"\t%lg",t->energy);
  fputc('\n',f);
  if (sw->seqdisp) disp_seq(f,t,sw); }


void disp_batch_tmrna(FILE *f, gene *t, csw *sw)
{ int tpe,*sb,*se;
  char pos[50];
  static char perm[3] = "NP";
  position(pos,t,sw);
  fprintf(f,"M%c\t%s\t",perm[(t->asst == 0)?0:1],pos);
  tpe = t->tpe;
  sb = t->eseq + ASTEM1_EXT + t->tps;
  se = t->eseq + ASTEM1_EXT + tpe + 1; 
  while (ltranslate(se) == '*') 
   { se += 3;
     tpe += 3; }
  fprintf(f,"[%d,%d]\t",t->tps+1,tpe+1);
  while (sb < se)  
   { fputc(ltranslate(sb),f);
     sb += 3; }
  if (sw->energydisp)
   fprintf(f,"\t%lg",t->energy);
  fputc('\n',f);
  if (sw->seqdisp) disp_seq(f,t,sw); }



int find_tstems(int *s, int ls, trna_loop hit[], int nh, csw *sw)
{ int i,r,c,tstem,tl,tloop,ithresh1;
  int *s1,*s2,*se,*ss,*si,*sb,*sc,*sf,*sl,*sx,*template;
  double e,ec,energy,penalty,thresh2;
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
  static int template_trna[6] = 
   { 0x0100, 0x0002, 0x2000, 0x0220, 0x0000, 0x0000 };
  static int template_tmrna[6] = 
   { 0x0100, 0x0002, 0x2220, 0x0220, 0x0000, 0x0000 };
  i = 0;
  template = (sw->tmrna)?template_tmrna:template_trna;
  ithresh1 = (int)sw->tcthresh;
  thresh2 = sw->tthresh;
  ss = s + sw->loffset;
  si = ss + 4 - 1;
  sl = s + ls - sw->roffset + 5 + 3;
  r = template[*si++];
  r = (r >> 4) + template[*si++];
  r = (r >> 4) + template[*si++];
  while (si < sl)
   { r = (r >> 4) + template[*si++];
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


int find_dstems(int *s, int ls,
                trna_dloop hit[], int nh,
		double thresh, int tfold)
{ int i,j,k,dstem,dloop,ige[7];
  int *s1,*s2,*sd,*si,*sb,*se,*sc,*sf,*sl,*sg1,*sg2;
  unsigned int q,r,c;
  static unsigned int GG[6] =
   { 0x00, 0x00, 0x11, 0x00, 0x00, 0x00 };
  static int qbp[6][6] =
   { { 0,0,0,1,0,0 },
     { 0,0,1,0,0,0 },
     { 0,1,0,1,0,0 },
     { 1,0,1,0,0,0 },
     { 0,0,0,0,0,0 },
     { 0,0,0,0,0,0 } };
  double genergy,energy,energy2,energyf,energyf6;
  static double bem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double fem[6][6] =
   { { -4.000,-4.000,-4.000, ATBOND, 0.000, 0.000 },
     { -4.000,-4.000, 3.000,-4.000, 0.000, 0.000 },
     { -4.000, 3.000,-4.000, 1.286, 0.000, 0.000 },
     {  ATBOND,-4.000, 1.286,-4.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double GC[6] = { 0.0,1.5,6.0,0.0,0.0,0.0 };
  static double G3[6] = { 0.0,6.0,12.0,12.0,0.0,0.0 };
  static double R[6] = { 6.0,0.0,6.0,0.0,0.0,0.0 };
  static double RH[6] = { 3.0,0.0,3.0,0.0,0.0,0.0 };
  static double AGT[6] = { 6.0,0.0,6.0,6.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,6.0,0.0,0.0 };
  static int goffb[13] = { 0,0,0,0,1,2,2,2,2,2,2,2,2 };
  static int goffe[13] = { 0,0,0,0,2,3,4,4,5,6,6,6,6 };
  i = 0;
  sc = s;
  energyf = fem[sc[5]][tfold];
  sl = s + ls;
  while (sc < sl)
   { energy2 = T[sc[-2]] + RH[*(sc-1)] + GC[*sc] + fem[sc[-2]][sc[4]];
     energyf6 = fem[sc[6]][tfold];
     for (dstem = 3; dstem <= 4; dstem++)
      { sd = sc + dstem;
        dloop = 3;
        se = sd + dloop;
        energy = energy2 + 6.0 + R[*(se-1)] + energyf;
        if (dstem == 3) 
         if (energyf < 0.0) energyf = energyf6;
        se += dstem;
	s1 = sc;
	s2 = se;
	sf = s1 + dstem;
	while (s1 < sf) energy += bem[*s1++][*--s2];
        if (energy >= thresh)
         { if (i >= nh) goto FL;
           hit[i].pos = sc;
	   hit[i].end = ++se;
           hit[i].loop = dloop;
           hit[i].stem = dstem;
           hit[i].energy = energy;
           i++; }
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
           energy = energy2 + genergy + R[*(se-1)] + energyf;
           se += dstem;
	   s1 = sc;
	   s2 = se;
	   sf = s1 + dstem;
	   while (s1 < sf) energy += bem[*s1++][*--s2];
           if (energy >= thresh)
            { if (i >= nh) goto FL;
	      hit[i].pos = sc;
	      hit[i].end = ++se;
              hit[i].loop = dloop;
              hit[i].stem = dstem;
	      hit[i].energy = energy;
              i++; }}}
      s1 = sc;
      s2 = sc + 16;
      sd = sc + 6;
      j = qbp[*s1][*--s2];
      while (++s1 < sd) j += qbp[*s1][*--s2];
      if (j >= 6)
       { energy = T[sc[-1]] + RH[*sc] + GC[*(sc+1)] + energyf6;
         energy += G[*++sd];
         energy += G[*++sd];
         energy += AGT[*++sd] + fem[sc[-1]][sc[4]];
         sd += 7;
         s1 = sc;
         s2 = sd;
         sf = s1 + 6;
         while (s1 < sf) energy += bem[*s1++][*--s2];
         if (energy >= thresh)
          { if (i >= nh) goto FL;
            hit[i].pos = sc;
            hit[i].end = ++sd;
            hit[i].loop = 4;
            hit[i].stem = 6;
            hit[i].energy = energy;
            i++; }}
      s1 = sc;
      s2 = sc + 18;
      sd = sc + 7;
      j = qbp[*s1][*--s2];
      while (++s1 < sd) j += qbp[*s1][*--s2];
      if (j >= 7)
       { energy = energy2 + fem[sc[7]][tfold];
         energy += G[*++sd];
         energy += G[*++sd];
         energy += AGT[*++sd];
         sd += 8;
         s1 = sc;
         s2 = sd;
         sf = s1 + 7;
         while (s1 < sf) energy += bem[*s1++][*--s2];
         if (energy >= thresh)
          { if (i >= nh) goto FL;
            hit[i].pos = sc;
            hit[i].end = ++sd;
            hit[i].loop = 4;
            hit[i].stem = 7;
            hit[i].energy = energy;
            i++; }}
     energyf = energyf6;
     sc++; }
  FN:
  return(i);
  FL: 
  fprintf(stderr,"Too many D-stem hits\n");
  return(i); }


int find_astem5(int *si, int *sl, int *astem3, int n3,
                trna_loop hit[], int nh, csw *sw)
{ int i,k;
  int *s1,*s2,*se;
  unsigned int r,itathresh;
  double tathresh,energy;
  static unsigned int template[6] = { 0,0,0,0,0,0 };
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
  itathresh = sw->itathresh;
  tathresh = sw->tathresh;
  i = 0;
  sl += n3;
  se = astem3 + n3 - 1;
  template[0] = A[*se];
  template[1] = C[*se];
  template[2] = G[*se];
  template[3] = T[*se];
  while (--se >= astem3)
   { template[0] = (template[0] << 4) + A[*se];
     template[1] = (template[1] << 4) + C[*se];
     template[2] = (template[2] << 4) + G[*se];
     template[3] = (template[3] << 4) + T[*se]; }
  r = template[*si++];
  k = 1;
  while (++k < n3) r = (r >> 4) + template[*si++];
  while (si < sl)
   { r = (r >> 4) + template[*si++];
     if ((r & 15) >= itathresh)
      { s1 = astem3;
        s2 = si;
        se = s1 + n3;
        energy = abem[*s1++][*--s2];
        while (s1  < se)
         energy += abem[*s1++][*--s2];
        if (energy >= tathresh)
         { if (i >= nh)
            { fprintf(stderr,"Too many astem5 hits\n");
              goto FN; }
           hit[i].pos = si - n3;
           hit[i].energy = energy;
           i++; }}}
  FN:
  return(i); }


int find_icstems(int *si, int *sl, trna_loop hit[], int nh)
{ int i,cloop;
  int *s1,*s2,*sb,*se,*sc;
  unsigned int r;
  double energy;
  static unsigned int template[6] = { 0,0,0,0,0,0 };
  static unsigned int A[6] = { 0,0,0,2,0,0 };
  static unsigned int C[6] = { 0,0,2,0,0,0 };
  static unsigned int G[6] = { 0,2,0,1,0,0 };
  static unsigned int T[6] = { 2,0,1,0,0,0 };
  static double cbem[6][6] =
   { { -1.072,-0.214,-1.072,2.0*ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 6.000,-1.072, 0.000, 0.000 },
     { -1.072, 6.000,-1.072, 3.400, 0.000, 0.000 },
     {  2.0*ATBOND,-1.072, 3.400,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  i = 0;
  sb = si;
  sc = si + 5;
  se = si + 4;
  template[0] = A[*se];
  template[1] = C[*se];
  template[2] = G[*se];
  template[3] = T[*se];
  while (--se >= si)
   { template[0] = (template[0] << 4) + A[*se];
     template[1] = (template[1] << 4) + C[*se];
     template[2] = (template[2] << 4) + G[*se];
     template[3] = (template[3] << 4) + T[*se]; }
  si += 11;
  se = sl - VARDIFF - 5;
  if (si < se) si = se;
  r = template[*si++];
  r = (r >> 4) + template[*si++];
  r = (r >> 4) + template[*si++];
  r = (r >> 4) + template[*si++];
  while (si < sl)
   { r = (r >> 4) + template[*si++];
     if ((r & 15) >= 5)
      { if (i >= nh)
         { fprintf(stderr,"Too many cstem hits\n");
           goto FN; }
        hit[i].pos = si;
        hit[i].stem = 5;
        hit[i].loop = (int)(si - sc - 5);
        s1 = sb;
	s2 = si;
	se = s1 + 5;
        hit[i].energy = cbem[*s1++][*--s2];
        while (s1  < se)
         hit[i].energy += cbem[*s1++][*--s2];
        i++; }}
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
  static unsigned int template[6] =
   { 0x10310000, 0x01000101, 0x00010030,
     0x02000100, 0x00000000, 0x00000000 };
  static int A[6] = { 0,1,1,1,1,1 };
  /* static int R[6] = { 0,1,0,1,1,1 }; */
  static int V[6] = { 0,0,0,1,1,1 };
  static int M[6] = { 0,0,1,1,1,1 };
  thresh = sw->tmrthresh;
  i = 0;
  sl = s + ls;
  r = template[*s++];
  r = (r >> 4) + template[*s++];
  r = (r >> 4) + template[*s++];
  r = (r >> 4) + template[*s++];
  r = (r >> 4) + template[*s++];
  r = (r >> 4) + template[*s++];
  r = (r >> 4) + template[*s++];
  if (sw->tmstrict)
    while (s < sl)
     { r = (r >> 4) + template[*s++];
       if ((c = (r & 0xF)) < thresh) continue;
       c -= (V[s[1]] + V[s[2]] + M[s[5]] + A[s[8]]);
       if (c < thresh) continue;
       if (i >= nh) goto FL;
       st = s - 2;
       si = st;
       sb = st + MINTAGDIST + 2;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (*si++ != 3) 
           si++;
          else
           if (*si == 0)
            { if (!(*++si & 5)) goto ST1; }
           else 
            if (*si == 2)
             { if (*++si == 0) goto ST1; } 
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
     { r = (r >> 4) + template[*s++];
       if ((c = (r & 0xF)) < thresh) continue;
       if (i >= nh) goto FL;
       st = s - 2;
       si = st + MINTAGDIST;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (*si++ != 3) 
           si++;
          else
           if (*si == 0)
            { if (!(*++si & 5)) goto ST2; }
           else 
            if (*si == 2)
             { if (*++si == 0) goto ST2; } 
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






gene *find_nearest_gene(data_set *d, gene *ts, int nt, gene *t, csw *sw)
{ int n,i,comp,status;
  long a,b,c,e,score,smax,thresh,psmax;
  static long proximity = 3*MINCTRNALEN/4;
  psmax = d->psmax;
  comp = t->comp;
  status = t->status;
  smax = -1;
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
           if (a <= e)
            if (b >= c)
             if (status == ts[i].status)
              if ((comp == ts[i].comp) || ((status == 0) && sw->mt))
               { score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
       	         if (score > smax)
	          { n = i;
	            smax = score; }}
           c -= psmax; 
           e -= psmax; }
        if (a <= e)
         if (b >= c)
          if (status == ts[i].status)
           if ((comp == ts[i].comp) || ((status == 0) && sw->mt))
            { score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
              if (score > smax)
               { n = i;
                 smax = score; }}}
     a -= psmax; 
     b -= psmax; }
  for (i = 0; i < nt; i++)
   { c = ts[i].start;
     e = ts[i].stop;
     if (e < c)
      { e += psmax;
        if (a <= e)
         if (b >= c)
          if (status == ts[i].status)
           if ((comp == ts[i].comp) || ((status == 0) && sw->mt))
            { score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
       	      if (score > smax)
	       { n = i;
	         smax = score; }}
        c -= psmax; 
        e -= psmax; }
     if (a <= e)
      if (b >= c)
       if (status == ts[i].status)
        if ((comp == ts[i].comp) || ((status == 0) && sw->mt))
         { score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?e-c:b-c);
           if (score > smax)
            { n = i;
              smax = score; }}}
  if (t->tps > 0)
   { if ((10*smax) > (9*thresh)) return(ts + n);
     else return(NULL); }
  if (++smax >= proximity) return(ts + n);
  return(NULL); }



void overlap(data_set *d, gene ts[], int sort[], int n, int it, csw *sw)
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




gene *find_slot(data_set *d, gene *t, gene *ts, int *nts, csw *sw)
{ char s1[80],s2[80],s3[80],s4[80];
  gene *tn,*tsr;
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
  tn = find_nearest_gene(d,ts,*nts,t,sw);
  if (tn)
   { if (t->energy <= tn->energy) return(NULL);
     copy(tn->name,t->name);
     if (sw->verbose)
      { fprintf(stderr,"%s %s ",name(t,s1,0,sw),position(s3,t,sw));
	if (sw->energydisp) fprintf(stderr,"(%lg) ",t->energy);
        fprintf(stderr,"replacing %s %s",name(tn,s2,1,sw),
		position(s4,tn,sw));
	if (sw->energydisp) fprintf(stderr," (%lg)",tn->energy);
	fprintf(stderr,"\n"); }}
  else
   { if (*nts >= sw->space)
      { sw->space += NT;
        tsr = (gene *)realloc((void *)ts,sw->space*sizeof(gene));
        if (tsr == NULL)
         { fprintf(stderr,"No more memory to store detected genes\n");
           fprintf(stderr,"Gene lost\n");
           sw->space -= NT; 
	   return(NULL); }
        ts = tsr; }
     copy3cr(d->seqname,t->name,79);
     tn = ts + (*nts);
     *nts = (*nts) + 1;
     if (sw->verbose)
      { fprintf(stderr,"%s at %s",name(t,s1,0,sw),position(s2,t,sw));
	if (sw->energydisp) fprintf(stderr," (%lg)",t->energy);
	fprintf(stderr,"\n"); }}
  return(tn); }




int tmopt(data_set *d,
          trna_loop *th, int tarm, double the,
          trna_loop *ahit, int nah,
	  gene *to, int nto,
	  int *seq, csw *sw)
{ int r,na,nr,nrh,ibase,flag,tmstrict;
  int *s,*v,*s1,*s2,*sa,*sb,*se,*sf,*ps,*tpos,pseq[MAXETRNALEN+1];
  static int gtemplate[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
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
  trna_loop rhit[NH];
  gene te,*tn;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,2,0,13,8,0,28,0,0,3,5,7,1,0.0,0 };
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
        r = gtemplate[*sb++];
        while (sb < sf)
         { r = (r >> 4) + gtemplate[*sb++];
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
   { te.nbase += ASTEM2_EXTD;
     te.start = (long)(te.ps - seq);
     tn = find_slot(d,&te,to,&nto,sw);
     if (tn)
      { te.intron = te.astem1 + te.dloop + te.cstem;
        te.asst = 0;
        base_copy3(te.ps,te.eseq,te.nbase + te.nintron);
        remove_intron(te.ps,pseq,te.nbase,
                      te.intron,te.nintron);
        base_copy3(pseq,te.seq,te.nbase);
        *tn = te; }}
  return(nto); }


int tmopt_perm(data_set *d,
          trna_loop *th, int tarm, double the,
          trna_loop *ahit, int nah,
	  gene *to, int nto,
	  int *seq, csw *sw)
{ int r,na,nr,nrh,ibase,flag,tmstrict;
  int *s,*v,*s1,*s2,*sa,*sb,*se,*sf,*ps,*apos,*tpos,pseq[MAXETRNALEN+1];
  static int gtemplate[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
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
  trna_loop rhit[NH];
  gene te,*tn;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,2,0,13,8,0,28,0,0,3,5,7,1,0.0,0 };
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
     r = gtemplate[*sb++];
     while (sb < sf)
      { r = (r >> 4) + gtemplate[*sb++];
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
                 te.astem2 + ASTEM2_EXT;
     base_copy3(te.ps-54,te.eseq,te.tpe+1+TMPTRAILER);
     te.nbase = te.astem1 + te.dloop + te.cstem;
     base_copy3(te.ps+te.asst,te.seq,te.nbase);
     base_copy3(te.ps,te.seq+te.nbase,te.intron);
     te.nbase += te.intron;
     te.nintron = te.tpe - te.nbase + 1 + TMPTRAILER;
     te.intron += 54;
     te.tps += 54;
     te.tpe += 54;
     te.asst += 54;
     tn = find_slot(d,&te,to,&nto,sw);
     if (tn) *tn = te; }
  return(nto); }



int tmioptimise(data_set *d, int *seq, int lseq, gene *to, int nto, csw *sw)
{ int i,j,intron,nt,nth,nd1,nd2,ndx,ndh,na,nah,nppah,nc,nch,tfold,tarm;
  int dstem,flag,mindist,maxdist,tmindist,tmaxdist,tmmindist,tmmaxdist;
  int tthresh,tmstrict,*pe,*s,pseq[2*MAXETRNALEN+1];
  int *se,*sc,*sb,*si,*tpos,*tend,*apos,*dpos,*tloopfold,*tmv,*cpos,*cend;
  unsigned int r;
  static unsigned int TT[6] =
   { 0x00, 0x00, 0x00, 0x11, 0x00, 0x00 };
  static unsigned int GG[6] =
   { 0x00, 0x00, 0x11, 0x00, 0x00, 0x00 };
  static int qbp[6][6] =
   { { 0,0,0,1,0,0 },
     { 0,0,1,0,0,0 },
     { 0,1,0,1,0,0 },
     { 1,0,1,0,0,0 },
     { 0,0,0,0,0,0 },
     { 0,0,0,0,0,0 } };
  static int yic[9] = { 1,0,0,0,0,0,0,0,0 };
  static int tic[9] = { 1,1,0,0,0,0,0,0,0 };
  static int a1ic[9] = { 1,1,1,0,0,0,0,0,0 };
  static int a2ic[9] = { 1,1,1,1,0,0,0,0,0 };
  static int a3ic[9] = { 1,1,1,1,1,0,0,0,0 };
  static int ric[9] = { 1,1,1,1,1,1,0,0,0 };
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
  double e,ec,f,he,the,thet,them,emax,energy,cenergy,denergy,ienergy;
  trna_loop thit[NTH],chit[NC],ahit[NA];
  trna_dloop dhit[ND];
  gene te,*tn;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,2,3,9,5,7,0,0,0,15,5,7,0,0.0,0 };
  if (sw->mt & (!sw->tmrna)) goto MT; 
  emax = sw->trnathresh;
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
  tthresh = sw->tthresh;
  tmstrict = sw->tmstrict;
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
     te.energy = emax;
     the = thit[nt].energy;
     nah = find_astem5(tpos-maxdist,tpos-mindist,tend,7,ahit,NA,sw);
     if (sw->tmrna)
      { thet = the - G[tpos[t.tstem]] - G[tpos[t.tstem+1]];
        if (tmstrict)
         { if (thet >= tthresh)
            nto = tmopt(d,thit+nt,tarm,thet,ahit,nah,to,nto,seq,sw); }
        else
         nto = tmopt(d,thit+nt,tarm,the,ahit,nah,to,nto,seq,sw); 
        nppah = find_astem5(tpos+MINPPASDIST,tpos+MAXPPASDIST,
                            tend,7,ahit+nah,NA-nah,sw);
        nto = tmopt_perm(d,thit+nt,tarm,the,ahit+nah,nppah,
                         to,nto,seq,sw);
        if (thet < tthresh) continue; 
        the = thet; } 
     if (!sw->trna) continue;
     na = -1;
     while (++na < nah)
      { apos = ahit[na].pos;
        if (apos < (tpos - tmaxdist)) continue;
        if (apos > (tpos - tmindist)) break;
        he = the + ahit[na].energy;
	ndh = find_dstems(apos+8,2,dhit,ND,sw->dthresh,tfold);
        nd1 = ndh;
        while (--nd1 >= 0)
	 { dstem = dhit[nd1].stem;
           dpos = dhit[nd1].pos;
           if ((int)(dpos - apos) < 9) dhit[nd1].energy -= 3.0;
	   if (*tloopfold == 2)
	    { sb = dpos + dstem + 2;
              sc = sb;
	      se = sb + t.dloop - 3;
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
                 ndx = nd2; }}
	   denergy += he;
           if (denergy < (te.energy - 45.0)) goto NX;
           nch = find_icstems(cpos,tmv,chit,NC);
	   nc = -1;
           while (++nc < nch)
            { energy = denergy + chit[nc].energy;
              if (energy < (te.energy - 15.0)) continue;
              cend = chit[nc].pos;
              t.var = (int)(tpos - cend);
              t.cloop = chit[nc].loop;
              t.cstem = chit[nc].stem;
              intron = 0;
              if (t.cloop < 9)
               { if (sw->minintronlen > 0) continue;
                 t.nintron = 0;
                 if (t.var > 17) energy += vloop_stability(cend,t.var);
		 sb = cpos + t.cstem;
                 energy += T[*(sb + 1)] + Y[*(sb)] + R[*(sb + 5)] - 
                           0.05*t.var - ((t.cloop == 7)?0.0:6.0); }
              else
               { t.nintron = t.cloop - 7;
                 if (t.nintron > sw->maxintronlen) continue;
	         if (t.nintron < sw->minintronlen) continue;
                 if (t.var > 17) energy += vloop_stability(cend,t.var);
                 if (energy < (te.energy - 5.0)) continue; 
                 t.cloop = 7;
                 sb = cpos + t.cstem;
                 se = sb + t.nintron;
                 cenergy = YP[*se] + T[*(se+1)] + RP[*(se+5)];
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
                       intron = j; }}
                  energy +=  cenergy - 10.0 - ilw*(t.nintron  + 1.1*t.var);
                  if (t.nintron >= 130)
                   { si = se + intron;
                     j = si[-1];
                     if (j != 2) 
                      { if (si[-2] != 0) energy -= 4.0;
                        if (j != 1)
                         if (j != 3)
                          energy -= 8.0; }}}
              dstem = dhit[ndx].stem;
              dpos = dhit[ndx].pos;
              if (dstem >= 6)
               { if (sb[2 + a1ic[intron]*t.nintron] != 3) continue; 
                 if (sb[3 + a2ic[intron]*t.nintron] != 1) continue;
                 if (sb[4 + a3ic[intron]*t.nintron] != 0) continue; 
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
	      if (energy >= te.energy)
               { flag = 1;
                 t.energy = energy;
                 t.dstem = dstem;
                 t.dloop = dhit[ndx].loop;
                 t.astem1 = (t.dstem < 6)?7:((t.tstem < 5)?9:8);
                 t.astem2 = t.astem1;
                 t.spacer = dpos - apos - 7;
	         t.ps = dpos - t.spacer - t.astem1;
	         j = (int)(cpos - t.ps) + t.cstem;
	         t.anticodon = j + 2;
                 te = t;
	         if (te.nintron > 0) te.intron = j + intron; }}
           NX:
	   nd2 = nd1;
	   while (--nd2 >= 0)
	    if (dhit[nd2].end == cpos) dhit[nd2].end = NULL; }}
     if (flag)
      { te.nbase = te.astem1 + te.spacer + 1 + 2*te.dstem +
                   te.dloop +  2*te.cstem + te.cloop +
                   te.var + tarm + te.astem2;
        if (te.astem1 == 7)
         if (qbp[te.ps[-1]][te.ps[te.nbase + te.nintron]])
         { te.ps--;
           te.nbase += 2;
           te.anticodon++;
           if (te.nintron > 0) te.intron++;
           te.astem1 = 8;
           te.astem2 = 8; }
	te.nbase += ASTEM2_EXTD;
        te.start = (long)(te.ps - seq);
        tn = find_slot(d,&te,to,&nto,sw);
	if (!tn) continue;
	if (te.nintron == 0)
          base_copy3(te.ps,te.seq,te.nbase + ASTEM2_EXTE);
        else
         { base_copy3(te.ps,te.eseq,te.nbase + te.nintron + ASTEM2_EXTE);
           remove_intron(te.ps,pseq,te.nbase+ASTEM2_EXTE,
                         te.intron,te.nintron);
           base_copy3(pseq,te.seq,te.nbase + ASTEM2_EXTE); }
	*tn = te; }}
  MT:
  return(nto); }



void init_gene(gene *t, int nstart, int nstop)
{ int i;
  for (i = nstart; i < nstop; i++)
   { t[i].energy = -1.0;
     t[i].status = -1;
     t[i].tps = 0;
     *(t[i].name) = '\0'; }}


void disp_freq_table(gene *t, int nt, csw *sw)
{ int i,j,k,m,*s,c[3],n[3],table[6][6][6];
  static int cgflip[4] = { 0,2,1,3 };
  FILE *f = sw->f;
  for (i = 0; i < 4; i++) 
   for (j = 0; j < 4; j++)
    for (k = 0; k < 4; k++)
     table[i][j][k] = 0;
  for (i = 0; i < nt; i++)
   if (t[i].energy >= 0.0)
    if (t[i].status == 0)
     { s = t[i].seq + t[i].anticodon;
       table[*s][s[1]][s[2]]++; }
  fprintf(f,"tRNA Anticodon Frequency\n");
  for (j = 0; j < 4; j++)
   { n[2] = cgflip[j];
     for (k = 0; k < 4; k++)
      { n[1] = cgflip[k];
        for (i = 0; i < 4; i++) 
         { n[0] = cgflip[i];
           fprintf(f,"%c%c%c",cpbase(n[0]),cpbase(n[1]),cpbase(n[2]));
           m = table[n[0]][n[1]][n[2]];
           if (m > 0)
            fprintf(f," %-4s %-5d",aa(n,sw),m);
           else
            fprintf(f," %-4s      ",aa(n,sw)); }
        fputc('\n',f); }}
  fputc('\n',f); 
  fprintf(f,"tRNA Codon Frequency\n");
  for (i = 0; i < 4; i++)
   { n[0] = 3 - cgflip[i];
     for (j = 0; j < 4; j++)
      { n[1] = 3 - cgflip[j];
        for (k = 0; k < 4; k++) 
         { n[2] = 3 - cgflip[k];
           fprintf(f,"%c%c%c",cpbase(n[0]),cpbase(n[1]),cpbase(n[2]));
           c[0] = 3 - n[2];
           c[1] = 3 - n[1];
           c[2] = 3 - n[0];
           m = table[c[0]][c[1]][c[2]];
           if (m > 0)
            fprintf(f," %-4s %-5d",aa(c,sw),m);
           else
            fprintf(f," %-4s      ",aa(c,sw)); }
        fputc('\n',f); }}
  fputc('\n',f); }

void disp_energy_stats(gene *t, int nt, csw *sw)
{ int i,n[NS],status,introns,nintron;
  static char status_name[NS][30] =
   { "tRNA genes","tmRNA genes","","","Overall" };
  FILE *f = sw->f;
  if ((sw->trna) && (sw->maxintronlen > 0))
   { introns = 1;
     nintron = 0; }
  else
   introns = 0;
  for (i = 0; i < NS; i++) n[i] = 0;
  for (i = 0; i < nt; i++)
   if (t[i].energy >= 0.0)
    { n[NS-1]++;
      status = t[i].status;
      n[status]++; 
      if (status == 0) if (introns) if (t[i].nintron > 0) nintron++; }
  fputc('\n',f);
  fputc('\n',f);
  if (sw->trna) 
   { sw->ngene[0] += n[0]; 
     if (n[0] > 3) disp_freq_table(t,nt,sw);
     if ((n[0] > 0) || ((sw->tmrna) && (n[1] > 0)))
      if (introns)
       { fprintf(f,"Number of tRNA genes with no introns = %d\n",
                 n[0]-nintron);
         fprintf(f,"Number of tRNA genes with C-loop introns = %d\n",
                 nintron); }
      else
       fprintf(f,"Number of %s = %d\n",status_name[0],n[0]); }
  if (sw->tmrna)
   { sw->ngene[1] += n[1]; 
     if ((n[1] > 0) || ((sw->trna) && (n[0] > 0)))
      fprintf(f,"Number of %s = %d\n",status_name[1],n[1]); }
  for (i = 2; i < NS-1; i++)
   if (n[i] > 0)
    { sw->ngene[i] += n[i]; 
      fprintf(f,"Number of %s = %d\n",status_name[i],n[i]); }
  fputc('\n',f);
  fputc('\n',f); }

int gene_sort(data_set *d, gene *ts, int nt, int sort[])
{ int i,n,j,k;
  long starti,startj,stopi,stopj;
  n = 0;
  for (i = 0; i < nt; i++)
   if (ts[i].energy >= 0.0)
    sort[n++] = i;
  i = -1;
  while (++i < (n-1))
   { j = i;
     while (++j < n)
      { starti = ts[sort[i]].start;
	startj = ts[sort[j]].start;
      if (starti > startj)
       { k = sort[i];
	 sort[i] = sort[j];
	 sort[j] = k; }
      else
       if (starti == startj)
	{ stopi = ts[sort[i]].stop;
          stopj = ts[sort[j]].stop;
	  if (stopi < starti) stopi += d->psmax;
	  if (stopj < startj) stopj += d->psmax;
          if (stopi > stopj)
           { k = sort[i];
	     sort[i] = sort[j];
	     sort[j] = k; }}}}
   return(n); }


void disp_gene_set(data_set *d, gene *ts, int nt, csw *sw)
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
  n = gene_sort(d,ts,nt,sort);
  j = sw->tmrna_struct[51];
  for (i = 52; i <= 57; i++) j += sw->tmrna_struct[i];
  if (j != ((sw->tmrna_struct[0] << 4) + 9)) return;
  if (sw->libflag < 2)
   { if (n > 0)
      for (j = 0; j < n;)
       { i = sort[j++];
         t = ts + i;
         switch(t->status)
          { case 0: init_matrix(m);
                    disp_gene(t,m,sw);
                    sprintf(s,"%d.",j);
                    xcopy(m,0,32,s,length(s));
                    disp_matrix(f,m,MATY); 
                    overlap(d,ts,sort,n,i,sw);
                    if (sw->seqdisp) disp_seq(f,t,sw);
                    if (t->nintron > 0) disp_intron(f,t); 
                    break; 
            case 1: if (sw->trnadomaindisp == 1)
                     { init_matrix(m);
                       disp_gene(t,m,sw);
                       sprintf(s,"%d.",j);
                       xcopy(m,0,32,s,length(s));
                       disp_matrix(f,m,MATY); }
                    else
                     { fprintf(f,"\n%d.\n",j);
                       disp_location(t,sw,"Location"); }
                    overlap(d,ts,sort,n,i,sw);
                    if (t->asst == 0) disp_tmrna_seq(f,t,sw);
                    else disp_tmrna_perm_seq(f,t,sw); 
                    break; }}
     else 
      if (*(d->seqname) != '\0')
       fprintf(f,"\nNothing found in %s\n\n\n",d->seqname); 
      else
       fprintf(f,"\nNothing found\n\n\n"); }
  disp_energy_stats(ts,nt,sw); 
  if (nt > NT) free((void *)sort); }


void batch_gene_set(data_set *d, gene *ts, int nt, csw *sw)
{ int i,j,n,vsort[NT],*sort;
  char m[MATX][MATY],s[20];
  gene *t;
  FILE *f = sw->f;
  if (nt <= NT)
   sort = vsort;
  else
  {  sort = (int *)malloc(nt*sizeof(int));
     if (sort == NULL)
      { fprintf(stderr,"Not enough memory to sort genes\n");
        exit(1); }}
  n = gene_sort(d,ts,nt,sort);
  j = sw->tmrna_struct[51];
  for (i = 52; i <= 57; i++) j += sw->tmrna_struct[i];
  if (j != ((sw->tmrna_struct[0] << 4) + 9)) return;
  if (sw->libflag < 2)
   { fprintf(f,"%d\n",n);
     for (j = 0; j < n; j++)
      { t = ts + sort[j];
        switch(t->status)
         { case 0: disp_batch_trna(f,t,sw);
                   break;
           case 1: disp_batch_tmrna(f,t,sw);
                   break;
           default: break; }}}
  if (nt > NT) free((void *)sort); }


	



void iopt_fastafile(data_set *d, csw *sw)
{ int i,nt,ns,flag,len;
  int *s,*sf,*se,*sc,*swrap;
  int seq[LSEQ+WRAP+1],cseq[LSEQ+1],wseq[WRAP+1];
  long gap,start,rewind,drewind,psmax,tmaxlen,vstart,vstop;
  gene *ts;
  FILE *f = sw->f;
  ts = (gene *)malloc(NT*sizeof(gene));
  if (ts == NULL)
   { fprintf(stderr,"Not enough memory available to store detected genes\n");
     exit(1); } 
  sw->space = NT;
  init_tmrna(f,sw);
  fprintf(f,"\nPlease reference the following paper if you use this\n");
  fprintf(f,"program as part of any published research.\n\n"); 
  fprintf(f,"Laslett, D. and Canback, B. (2004) ARAGORN, a\n");
  fprintf(f,"program for the detection of transfer RNA and\n"); 
  fprintf(f,"transfer-messenger RNA genes in nucleotide sequences.\n");
  fprintf(f,"Nucleic Acids Research, 32;11-16.\n\n"); 
  fprintf(f,"\nSearching for ");
  if (sw->trna)
   { fprintf(f,"tRNA genes");
     if (sw->maxintronlen > 0) fprintf(f," with introns in anticodon loop");
     else fprintf(f," with no introns");
     if (sw->tmrna) fprintf(f," and tmRNA genes");
     fputc('\n',f);
     if (sw->maxintronlen > 0)
       fprintf(f,"Intron length from %d to %d bases\n",
            sw->minintronlen,sw->maxintronlen); 
     fprintf(f,"Using standard genetic code for anticodon");
     if (sw->tmrna) fprintf(f," and peptide tag");
     fprintf(f," prediction\n"); }
  else
  { fprintf(f,"tmRNA genes\n");
    fprintf(f,"Using standard genetic code for peptide tag prediction\n"); }
  if (sw->linear)
   fprintf(f,"Assuming linear topology, search will not wrap around ends\n");
  else
   fprintf(f,"Assuming circular topology, search wraps around ends\n");
  if (sw->both)
   fprintf(f,"Searching both strands\n");
  else
   fprintf(f,"Searching single strand only\n");
  fputc('\n',f);
  fputc('\n',f);
  rewind = MAXTAGDIST + 20;
  if (sw->trna) 
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
  ns = 0;
  d->nextseq = 0L;
  while (d->nextseq >= 0L)
   { d->seqstart = d->nextseq;
     if (!seq_init(d)) break;
     psmax = d->psmax;
     if (sw->verbose) 
      { fprintf(stderr,"%s\n",d->seqname);
        fprintf(stderr,"%ld nucleotides in sequence\n",psmax);
        fprintf(stderr,"Mean G+C content = %2.1lf%%\n",100.0*d->gc); }
     fprintf(f,"%s\n",d->seqname);
     fprintf(f,"%ld nucleotides in sequence\n",psmax);
     fprintf(f,"Mean G+C content = %2.1lf%%\n",100.0*d->gc);
     init_gene(ts,0,NT);
     nt = 0;
     flag = 0;
     start = 1L;
     se = seq;
     if (sw->linear)
      { for (i = 0; i < rewind; i++) *se++ = NOBASE;
        start -= rewind; } 
     else
      { if (psmax <= rewind)
	 { gap = rewind - psmax;
           sc = se + gap;
	   while (se < sc) *se++ = NOBASE;
           swrap = wseq;
           sc = se + psmax; 
           while (se < sc) 
            { *se = move_forward(d);
              *swrap++ = *se++; }
	   sc = swrap + gap;
	   while (swrap < sc) *swrap++ = NOBASE; 
           start -= gap; }
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
	    { sc = wseq + ((psmax <= rewind)?rewind:drewind);
              swrap = wseq;
              while (swrap < sc) *se++ = *swrap++; }
           flag = 1;
           break; }}
     len = (int)(se - seq);
     if (sw->verbose)
      { vstart = sq(start + rewind);
	vstop = sq(start + len - rewind - 1);
	if (vstop < vstart)
	 { fprintf(stderr,"Searching from %ld to %ld\n",vstart,psmax);
           fprintf(stderr,"Searching from 1 to %ld\n",vstop); }
	else
         fprintf(stderr,"Searching from %ld to %ld\n",vstart,vstop); }
     sw->start = start;
     sw->comp = 0;
     nt = tmioptimise(d,seq,len,ts,nt,sw);
     if (sw->both)
      { sense_switch(seq,cseq,len);
        sw->start = start+len;
        sw->comp = 1;
        nt = tmioptimise(d,cseq,len,ts,nt,sw); }
     if (!flag)
      {	s = seq;
	sf = se - drewind;
	se = seq + drewind;
	while (s < se) *s++ = *sf++;
        start += len - drewind;
        goto NX; }
     disp_gene_set(d,ts,nt,sw);
     ns++; }
  if (ns > 1)
   { fprintf(f,"\n\n%d sequences searched\n",ns);
     if (sw->trna) fprintf(f,"Total tRNA genes = %d\n",sw->ngene[0]);
     if (sw->tmrna) fprintf(f,"Total tmRNA genes = %d\n",sw->ngene[1]); }
  if (sw->verbose) fprintf(stderr,"Search Finished\n\n"); 
  free((void *)ts); }


void bopt_fastafile(data_set *d, csw *sw)
{ int i,nt,ns,flag,len;
  int *s,*sf,*se,*sc,*swrap;
  int seq[LSEQ+WRAP+1],cseq[LSEQ+1],wseq[WRAP+1];
  long gap,start,rewind,drewind,psmax,tmaxlen,vstart,vstop;
  gene *ts;
  FILE *f = sw->f;
  ts = (gene *)malloc(NT*sizeof(gene));
  if (ts == NULL)
   { fprintf(stderr,"Not enough memory available to store detected genes\n");
     exit(1); } 
  sw->space = NT;
  rewind = MAXTAGDIST + 20;
  if (sw->trna) 
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
  ns = 0;
  d->nextseq = 0L;
  while (d->nextseq >= 0L)
   { d->seqstart = d->nextseq;
     if (!seq_init(d)) break;
     psmax = d->psmax;
     fprintf(f,">%s\t",d->seqname);
     init_gene(ts,0,NT);
     nt = 0;
     flag = 0;
     start = 1L;
     se = seq;
     if (sw->linear)
      { for (i = 0; i < rewind; i++) *se++ = NOBASE;
        start -= rewind; } 
     else
      { if (psmax <= rewind)
	 { gap = rewind - psmax;
           sc = se + gap;
	   while (se < sc) *se++ = NOBASE;
           swrap = wseq;
           sc = se + psmax; 
           while (se < sc) 
            { *se = move_forward(d);
              *swrap++ = *se++; }
	   sc = swrap + gap;
	   while (swrap < sc) *swrap++ = NOBASE; 
           start -= gap; }
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
	    { sc = wseq + ((psmax <= rewind)?rewind:drewind);
              swrap = wseq;
              while (swrap < sc) *se++ = *swrap++; }
           flag = 1;
           break; }}
     len = (int)(se - seq);
     sw->start = start;
     sw->comp = 0;
     nt = tmioptimise(d,seq,len,ts,nt,sw);
     if (sw->both)
      { sense_switch(seq,cseq,len);
        sw->start = start+len;
        sw->comp = 1;
        nt = tmioptimise(d,cseq,len,ts,nt,sw); }
     if (!flag)
      {	s = seq;
	sf = se - drewind;
	se = seq + drewind;
	while (s < se) *s++ = *sf++;
        start += len - drewind;
        goto NX; }
     batch_gene_set(d,ts,nt,sw);
     ns++; }
  if (ns > 1)
   { fprintf(f,">end \t%d sequences",ns);
     if (sw->trna) fprintf(f," %d tRNA genes",sw->ngene[0]);
     if (sw->tmrna) fprintf(f," %d tmRNA genes",sw->ngene[1]);
     fputc('\n',f); }
  free((void *)ts); }


void aragorn_help_menu()
{ printf("\n");
  printf("---------------------------\n");
  printf("ARAGORN v1.1 Dean Laslett\n");
  printf("---------------------------\n");
  printf("\n");
  printf("Please reference the following paper if you use this\n");
  printf("program as part of any published research.\n\n"); 
  printf("Laslett, D. and Canback, B. (2004) ARAGORN, a\n");
  printf("program for the detection of transfer RNA and transfer-messenger\n"); 
  printf("RNA genes in nucleotide sequences\n");
  printf("Nucleic Acids Research, 32;11-16\n"); 
  printf("ARAGORN detects tRNA and tmRNA genes.\n");
  printf("\n");
  printf("Usage:  aragorn -v -s -d -c -l -i<min>,<max> -t -m -o <outfile>");
  printf(" <filename>\n\n");
  printf("<filename> is assumed to contain one or more sequences\n");
  printf("in FASTA format. Results of the search are printed to\n");
  printf("STDOUT. All switches are optional and case-insensitive.\n");
  printf("Unless -i is specified, tRNA genes containing introns\n");
  printf("are not detected. \n");
  printf("\n");
  printf("    -m            Search for tmRNA genes only.\n");
  printf("    -t            Search for tRNA genes only.\n");
  printf("                  By default, both are detected.\n");
  printf("    -i            Search for tRNA genes with introns in\n");
  printf("                  anticodon loop with maximum length %d\n",
         MAXINTRONLEN);
  printf("                  bases. Minimum intron length is 0 bases.\n");
  printf("                  Ignored if -m is specified.\n");
  printf("    -i<max>       Search for tRNA genes with introns in\n");
  printf("                  anticodon loop with maximum length <max>\n");
  printf("                  bases. Minimum intron length is 0 bases.\n");
  printf("                  Ignored if -m is specified.\n");
  printf("    -i<min>,<max> Search for tRNA genes with introns in\n");
  printf("                  anticodon loop with maximum length <max>\n");
  printf("                  bases, and minimum length <min> bases.\n");
  printf("                  Ignored if -m is specified.\n");
  printf("    -c            Assume that each sequence has a circular\n");
  printf("                  topology. Search wraps around each end.\n");
  printf("                  Default setting.\n");
  printf("    -l            Assume that each sequence has a linear\n");
  printf("                  topology. Search does not wrap.\n");
  printf("    -d            Double. Search both strands of each\n");
  printf("                  sequence. Default setting.\n");
  printf("    -s            Single. Do not search the complementary\n");
  printf("                  strand of each sequence.\n");
  printf("    -seq          Print out primary sequence.\n");
  printf("    -v            Verbose. Prints out search progress\n");
  printf("                  to STDERR.\n");
  printf("    -a            Print out tRNA domain for tmRNA genes.\n");
  printf("    -o <outfile>  print output into <outfile>. If <outfile>\n");
  printf("                  exists, it is overwritten.\n");
  printf("                  By default, output goes to STDOUT.\n");
  printf("    -w            Print out genes in batch mode.\n");
  printf("                  For tRNA genes, output is in the form:\n");
  printf("                  sequence name<tab>Ng<ret>\n");
  printf("                  [locus 1]<tab>T or TI<tab><Apos><tab><species>(<ACT>)\n");
  printf("                  <tab>intron position<tab>intron length<ret>\n");
  printf("                            .          \n");
  printf("                  [locus Ng]<tab>T or TI<tab><Apos><tab><species>(<ACT>)\n");
  printf("                  <tab>intron position<tab>intron length<ret>\n");
  printf("                  Ng is the number of genes found\n");
  printf("                  T means the tRNA gene has no C-loop intron\n");
  printf("                  TI means the tRNA gene has a C-loop intron\n");
  printf("                  <Apos> is the tRNA anticodon relative position\n");
  printf("                  <species> is the tRNA iso-acceptor species\n");
  printf("                  <ACT> is the tRNA anticodon base triplet\n");
  printf("                  For tmRNA genes, output is in the form:\n");
  printf("                  sequence name<tab>Ng<ret>\n");
  printf("                  [locus 1]<tab>MN or MP<tab>[tag offset]\n");
  printf("                  <tab>tag peptide<ret>\n");
  printf("                            .          \n");
  printf("                  [locus Ng]<tab>MN or MP<tab>[tag offset]\n");
  printf("                  <tab>tag peptide<ret>\n");
  printf("                  Ng is the number of tmRNA genes found\n");
  printf("                  MN means the tmRNA gene is not permuted\n");
  printf("                  MP means the tmRNA gene is permuted\n");
  printf("\n\n"); }

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
     default: break; }
  exit(0); }


int main(int z, char *v[])
{ int i,j,l,filecounter;
  char c,*s;
  data_set d;
  static csw sw =
   { "",NULL,0,1,1,0,0,0,5,5,1,0,0,0,1,0,0,0,132.0,325.0,10.85,4.0,
     29.0,26.0,7.5,14.0,10.0,25.0,1,9,8,0,0,-7.0,-4.0,-10.0,-7.9,-7.9,
     {0,0,0,0,0},0,0,0,0,0,
     { 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 10,
       65, 82, 65, 71, 79, 82, 78, 32,
       118, 49, 46, 49, 32, 32, 32, 68, 101, 97,
       110, 32, 76, 97, 115, 108, 101, 116, 116, 10,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 10, TERM }};
  sw.f = stdout;
  filecounter = 0;
  i = 0;
  while (++i < z)
   if (*(v[i]) == '-')
    { s = v[i] + 1;
      c = upcasec(*s);
      switch(c)
       { case  'E': sw.energydisp = 1;
                    break;
         case  'A': sw.trnadomaindisp = 1;
                    break;
         case  'W': sw.batch = 1;
                    break;
         case  'V': sw.verbose = 1;
                    break;
         case  'S': if (upcasec(s[1]) == 'E')
                     sw.seqdisp = 1;
                    else
                     sw.both = 0;
                    break;
         case  'D': sw.both = 1;
                    break;
         case  'L': sw.linear = 1;
                    break;
         case  'C': sw.linear = 0;
                    break;
         case  '1': sw.minintronlen = 10;
                    break;
         case  'I': if (length(v[i]) > 2)
                     { s = iconvert(v[i] + 2,&l);
		       if (*s == ',')
		        { sw.minintronlen = l;
			  iconvert(s+1,&sw.maxintronlen); }
		       else
		        sw.maxintronlen = l;
		       if (sw.maxintronlen > (LSEQ - MAXTRNALEN))
		         sw.maxintronlen = (LSEQ - MAXTRNALEN);
		       if (sw.maxintronlen > MAXINTRONLEN)
		        sw.maxintronlen = MAXINTRONLEN;
		       if ((sw.minintronlen < 0) ||
		           (sw.maxintronlen < sw.minintronlen))
			error_report(1,v[i]); }
		    else
		     sw.maxintronlen = MAXINTRONLEN;
                    break;
         case  'T': sw.tmrna = 0;
                    if (length(v[i]) > 2)
                     { s = dconvert(v[i] + 2,&sw.trnathresh);
		       if (*s == ',')
         	        dconvert(s+1,&sw.tthresh); }
                    break;
         case  'M': sw.trna = 0;
                    if (length(v[i]) > 2)
                     dconvert(v[i] + 2,&sw.tmrnathresh);
                    break;
         case  'R': sw.tmstrict = 0;
                    break;
         case  'Q': sw.showconfig = 0;
                    break;
         case  'H': aragorn_help_menu();
                    exit(0);
         case  'O': if (length(v[i]) > 2)
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
     error_report(0,v[i]);
  if (filecounter < 1)
   error_report(3,NULL);
  if (sw.libflag) fprintf(sw.f,"tRNA Library\n");
  if (sw.batch) bopt_fastafile(&d,&sw);
  else iopt_fastafile(&d,&sw);
  fclose(d.f);
  if (!sw.batch && sw.showconfig)
   { fprintf(sw.f,"Configuration: ");
     i = -1;
     while (++i < z) fprintf(sw.f,"%s ",v[i]);
     fputc('\n',sw.f); }
  if (sw.f != stdout) fclose(sw.f);
  return(0); }


