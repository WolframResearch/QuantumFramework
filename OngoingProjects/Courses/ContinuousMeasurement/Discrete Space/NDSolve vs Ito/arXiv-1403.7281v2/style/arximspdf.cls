%% This is file `arximspdf.cls'
%%
%% LaTeX 2e class file for the processing of LaTeX2e files
%% of the following IMS journals:
%%
%%   The Annals of Probability
%%   The Annals of Applied Probability
%%   The Annals of Statistics
%%
%%   Institute of Mathematical Statistics, U.S.A.
%%
%% Macros written by Vytas Statulevicius, VTeX, Lithuania
%% for Institute of Mathematical Statistics, U.S.A.
%% Please submit bugs or your comments to vytas@vtex.lt
%%
%% These class file is a simplified version of style files,
%% used in a professional typesetting of the above mentioned journals.
%% The main purpose of this class is to let users compile files
%% downloaded from the arXiv 
%%
%% The following macro packages are required:
%%   textcomp.sty
%%   amssymb.sty
%%   amsfonts.sty
%%   bm.sty
%%   texsort.sty
%%   keyval.sty
%%   natbib.sty
%%
%% You are free to use this class file as you see fit, provided 
%% that you do not make changes to the file. 
%% If you DO make changes, you are required to rename this file.
%%
%% It may be distributed under the terms of the LaTeX Project Public
%% License, as described in lppl.txt in the base LaTeX distribution.
%% Either version 1.0 or, at your option, any later version.
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
%%
%% Bug fixes and changes:
%% 2006.03.25 - author's addresses added to bookmarks
%% 2006.05.02 - changes in \ead \printead
%% 2006.05.15 - fixed bug in \printead
%% 2006.05.30 - fixed bug in \ead
%% 2006.08.03 - commands \referstodoi{} and \relateddoi{}{} implemented
%% 2007.08.22 - support for The  Annals of Applied Statistics (aoas) added
%% 2008.03.06 - supplement environment introduced
%% 2008.09.17 - doi in supplement changed to DOI
%% 2008.10.20 - msnbibl macros added
%% 2008.12.11 - \tableofcontents added
%% 2009.12.10 - msnbibl macros fixed
%% 2009.12.16 - \betal added
%% 2010.10.05 - in authoryear case links are provided only from year (solves layout problem)
%% 2011.03.08 - \sfilename introduced
%% 2012.01.12 - bibliography markup commands added
%% 2012.02.07 - bibliography markup commands added
%% 2012.05.29 - vtexurl macros added; \bsuffix added (ES)
%% 2012.06.20 - paragraph indentation after end of theorem (ES)
%% 2012.07.17 - auto keywords afterdot (ES)
%% 2012.10.17 - \blocation added (VS)
%% 2013.02.22 - binbook title font changed (ES)
%% 2013.07.01 - thanksmark marks changed in affiliation for AOAS (ES)
%% 2014.05.15 - \relateddoi{}{} and its \doi{} changed to \relateddois{}{} and \relateddoi{}{} respectively
%% 2014.11.27 - (end)OrigBibText (ES)
%% 2014.12.23 - thanksmark marks changed for the rest juornals (ES)
%
% General class file for the journals:
%
% TeX programming: Vytas Statulevicius, VTeX, Lithuania, vytas@vtex.lt
% Requires Latex2e, ver.2000.06


\def\fmt@version{2014/12/23}

\NeedsTeXFormat{LaTeX2e}[2000/06/01]
\ProvidesClass{arximspdf}
              [\fmt@version, IMS]

% Definitions:
\newcommand\bibname{References}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{}
\newcommand\acknowledgementsname{Acknowledgements}

\newcommand\@ptsize{}
\newif\if@restonecol \@restonecolfalse
\newif\if@openright  \@openrightfalse

% General Options
\DeclareOption{10pt}{}
\DeclareOption{11pt}{}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{openright}{\@openrighttrue}
\DeclareOption{openany}{\@openrightfalse}
\DeclareOption{onecolumn}{\@twocolumnfalse\@restonecoltrue}
\DeclareOption{twocolumn}{\@twocolumntrue}

\newif\if@noMRlinks \@noMRlinksfalse
\DeclareOption{noMRlinks}{\@noMRlinkstrue}

% naudojame MathSciNet bibl:

\newif\if@MSNbibl \@MSNbiblfalse

\DeclareOption{MSNbibl}{\@MSNbibltrue}


% Put "." after inline section headings:
\newif\if@autosecdot 
\DeclareOption{autosecdot}  {\@autosecdottrue}
\DeclareOption{noautosecdot}{\AtBeginDocument{\@autosecdotfalse}}

\DeclareOption{aap}{
  \def\@journal{The Annals of Applied Probability}
  \def\@issn{1050-5164}
  \def\journal@url{http://www.imstat.org/aap/}
}

\DeclareOption{aop}{%
   \def\@journal{The Annals of Probability}
   \def\@issn{0091-1798}
   \def\journal@url{http://www.imstat.org/aop/}%
}

\DeclareOption{aos}{
    \def\@journal{The Annals of Statistics}
    \def\@issn{0090-5364}
    \def\journal@url{http://www.imstat.org/aos/}%
}


\DeclareOption{aoas}{%
   \def\@journal{The  Annals of Applied Statistics}
   \def\@issn{1932-6157}
   \def\@jid{AOAS}
   \def\journal@url{http://www.imstat.org/aoas/}
}


% Information about journals
\def\@issn{0000-0000}
\def\@journal{Undefined}
\def\journal@url{http://www.imstat.org}
\def\owner@url{http://www.imstat.org}
\def\msc@url{http://www.ams.org/msc/}

% bibliography
\newif\if@nameyear \@nameyearfalse 
\DeclareOption{nameyear}{\global\@nameyeartrue}
\DeclareOption{number}{\global\@nameyearfalse}

% vtexflt.sty
\def\float@numbering{0}
\DeclareOption{normalfloat}{\def\float@numbering{0}}
\DeclareOption{secfloat}   {\def\float@numbering{1}}
\DeclareOption{chapfloat}  {\def\float@numbering{2}}
\DeclareOption{rotating}{\AtEndOfClass{\RequirePackage[dvips,counterclockwise,figuresright]{rotating}}}%

% vtexmth.sty
\def\eqn@numbering{0}
\DeclareOption{normaleqn}{\def\eqn@numbering{0}}
\DeclareOption{seceqn}   {\def\eqn@numbering{1}}
\DeclareOption{chapeqn}  {\def\eqn@numbering{2}}

% vtexthm.sty
\def\thm@numbering{0}
\DeclareOption{normalthm}{\def\thm@numbering{0}}
\DeclareOption{secthm}   {\def\thm@numbering{1}}
\DeclareOption{chapthm}  {\def\thm@numbering{2}}
\newif\if@xx@theorem \@xx@theoremfalse 
\DeclareOption{xxtheorem}{\global\@xx@theoremtrue} 

% hyperref.sty
\newif\if@loadhyperref
\DeclareOption{nohyperref}{\global\@loadhyperreffalse}
\@loadhyperreftrue

\PassOptionsToPackage{unknownkeysallowed}{keyval}

\ExecuteOptions{leqno,twoside,final,number,autosecdot}
\ProcessOptions


%************************* FONTS

\typeout{Eleven point}
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xipt{13\p@ plus .5\p@ minus .5\p@}%
   \abovedisplayskip 7\p@ \@plus3\p@ \@minus3\p@%
   \belowdisplayskip\abovedisplayskip
   \abovedisplayshortskip\abovedisplayskip%
   \belowdisplayshortskip\abovedisplayskip}

\newcommand\small{%
   \@setfontsize\small\@xpt{12pt plus .3pt minus .3pt}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus2\p@%
   \belowdisplayskip \abovedisplayskip%
   \abovedisplayshortskip\abovedisplayskip%
   \belowdisplayshortskip\abovedisplayskip}

\newcommand\footnotesize{%
   \@setfontsize\small\@ixpt{11pt plus .2\p@ minus .2\p@}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus2\p@%
   \belowdisplayskip \abovedisplayskip%
   \abovedisplayshortskip\abovedisplayskip%
   \belowdisplayshortskip\abovedisplayskip}%

\newcommand\scriptsize{\@setfontsize\scriptsize\@viiipt\@xpt}
\newcommand\tiny{\@setfontsize\tiny\@vipt\@viipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{16}}
\newcommand\LARGE{\@setfontsize\LARGE\@xvipt{18}}
\newcommand\huge{\@setfontsize\huge\@xviiipt{20}}
\newcommand\Huge{\@setfontsize\Huge\@xxpt{24}}

% Skliaustu parinkimo pritaikymas 11pt
\def\big#1{{\hbox{$\left#1\vbox to9.31\p@{}\right.\n@space$}}}
\def\Big#1{{\hbox{$\left#1\vbox to12.59\p@{}\right.\n@space$}}}
\def\bigg#1{{\hbox{$\left#1\vbox to15.88\p@{}\right.\n@space$}}}
\def\Bigg#1{{\hbox{$\left#1\vbox to19.16\p@{}\right.\n@space$}}}

% Customization of fonts
\renewcommand\sldefault{it}
\renewcommand\bfdefault{b}
\let\slshape\itshape

% Compatibility with old latex:
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\let\sl\it
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}

% ********************* DIMENSIONS:
% TEXT DIMENSIONS
\setlength\parindent    {12\p@}
\def\sv@parindent       {12\p@}
\setlength\textwidth    {30pc}
\setlength\textheight   {550\p@} %46 eil.
\setlength\columnsep    {10\p@}
\setlength\columnwidth  {270\p@}
\setlength\columnseprule{0\p@}

% HEADS:
\setlength\headheight{12\p@}
\setlength\headsep   {14\p@}
\setlength\topskip   {10\p@}
\setlength\footskip  {14\p@}
\setlength\maxdepth  {.5\topskip}

% PAPER AND TRIM SIZE
\setlength\paperwidth  {8.5in}
\@settopoint\paperwidth
\setlength\paperheight {11in}
\@settopoint\paperheight

% CENTER ON PAGE
    \setlength\@tempdima        {\the\paperwidth}
    \addtolength\@tempdima      {-\textwidth}
    \setlength\oddsidemargin    {.5\@tempdima}
    \addtolength\oddsidemargin  {-1in}
    \evensidemargin=\oddsidemargin
    \@settopoint\evensidemargin
    \@settopoint\oddsidemargin

    \setlength\@tempdima        {\the\paperheight}
    \addtolength\@tempdima      {-\textheight}
    \setlength\topmargin    {.5\@tempdima}
    \addtolength\topmargin  {-1in}

\if@twocolumn
 \setlength\marginparsep {10\p@}
\else
  \setlength\marginparsep{7\p@}
\fi
\setlength\marginparpush {5\p@}

% TEXT PARAMETERS
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@}

% BREAKS
\setlength\smallskipamount{6.5\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount  {13\p@ \@plus 3\p@ \@minus 3\p@}
\setlength\bigskipamount  {19.5\p@ \@plus 3\p@ \@minus 3\p@}

% PAGE-BREAKING PENALTIES
\clubpenalty=10000
\widowpenalty=10000
\displaywidowpenalty=50
\predisplaypenalty=10000   % Breaking before a math display.
\hfuzz=1pt            % do not be to critical about boxes

% FOOTNOTES
\setlength\footnotesep   {7\p@}
\setlength{\skip\footins}{18\p@ \@plus 4\p@ \@minus 2\p@}
\skip\@mpfootins = \skip\footins

% FRAMED BOXES
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}

%************************* LOADING PACKAGES:

\normalfont\normalsize

\RequirePackage[psamsfonts]{amssymb}
\newcommand\hmmax{0}
\newcommand\bmmax{1}
\RequirePackage{bm}

% FRONT MATTER MACROS

%%%% start of vtexfm.sty [2003/03/13 VTEX front matter general macros (VS)]
%%%%%%%%%%%%%% Options
\newif\if@journal %\DeclareOption{journal}{\@journaltrue}

\@journaltrue

\RequirePackage{texsort}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setattribute, getattribute, do@option@list

\def\setattribute{\@ifnextchar[\@setattribute{\@setattribute[]}}
\def\@setattribute[#1]#2#3#4{\expandafter\gdef\csname #2@#3\endcsname{#4}}
\def\getattribute#1#2{\csname #1@#2\endcsname}
\def\sep@key@value#1=#2/?/#3{\setattribute{#3}{#1}{#2}}
\def\do@option@list#1#2{%
  \@for\curr@option:={#2}\do{%
    \expandafter\sep@key@value\curr@option/?/{#1}\relax
  }%
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% firstpage, lastpage, pagerange, pagenumbering

% Numeruojame nuo 1:
\def\firstpage#1{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \gdef\@firstpage{#1}\gdef\@lastpage{#1}%
%  \global\c@page=#1 
  \ignorespaces\fi}

\def\lastpage#1{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \gdef\@lastpage{#1}\ignorespaces\fi}

\def\@firstpage{1}
\def\@lastpage{0}

\def\pagerange@sep{--}
\def\set@pagerange{%
    \ifx\@firstpage\@lastpage%
       \gdef\@pagerange{\@firstpage}%
    \else%
       \gdef\@pagerange{\@firstpage\pagerange@sep\@lastpage}%
    \fi}
\def\@pagerange{1\pagerange@sep 0}

\def\pagenumbering#1{%
  \def\@tempa{\csname none\endcsname}%
  \def\@tempb{\csname #1\endcsname}%
  \ifx\@tempa\@tempb
    \global\let\thepage\relax
    \gdef\firstpage##1{\gdef\@firstpage{}\dbl@error{Su \string\pagenumbering{none} \string\firstpage\space nenaudojamas!}}%
    \gdef\lastpage##1{\gdef\@lastpage{}\dbl@error{Su \string\pagenumbering{none} \string\lastpage\space nenaudojamas!}}%
  \else
    \gdef\thepage{\csname @#1\endcsname\c@page}%
  \fi
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% journal, volume, issue, pubyear

\if@journal

% JOURNAL
\def\journal#1{\gdef\@journal{#1}}     
\@ifundefined{@journal}{\gdef\@journal{Journal not defined}}{}

% VOLUME
\setattribute{volume} {multisep}{--} %1,2,3 -> 1-3
\setattribute{volume} {twosep}  {--} %1,2   -> 1-2

\setattribute{volume} {multilabel} {Vols.} %1,2,3 -> Vols.
\setattribute{volume} {twolabel}   {Vols.} %1,2   -> Vols.
\setattribute{volume} {label}      {Vol.}  %1     -> Vol.

\def\volume#1{%
  \initarray{@tempa}{#1}\@tempcnta\getarraylenght{@tempa}\relax%
  \ifnum\@tempcnta>1\global\let\volume@label\volume@twolabel\fi%
  \ifnum\@tempcnta>2\global\let\volume@label\volume@multilabel\fi%
  \compresslist{#1}\to{multivol@list}\sep{\volume@twosep}{\volume@multisep}%
  \expandafter\gdef\expandafter\@volume\expandafter{\multivol@list}}

\def\@volume{0}

%% ISSUE
\setattribute{issue} {multisep} {--} %1,2,3 -> 1-3
\setattribute{issue} {twosep}   {--} %1,2   -> 1-2

\setattribute{issue} {multilabel} {Nos.} %1,2,3 -> Nos.
\setattribute{issue} {twolabel}   {Nos.} %1,2   -> Nos.
\setattribute{issue} {label}      {No.}  %1     -> No.

\def\issue#1{%
  \initarray{@tempa}{#1}\@tempcnta\getarraylenght{@tempa}\relax%
  \ifnum\@tempcnta=1 \gdef\@issue{#1}%
  \else%
    \ifnum\@tempcnta>1\global\let\issue@label\issue@twolabel\fi%
    \ifnum\@tempcnta>2\global\let\issue@label\issue@multilabel\fi%
    \compresslist{#1}\to{multiissue@list}\sep{\issue@twosep}{\issue@multisep}%
    \expandafter\gdef\expandafter\@issue\expandafter{\multiissue@list}%
  \fi%
  }



\def\@issue{0}

% PUBYEAR
% \@pubyear 2001
% \copyyear   01

\newcount\@pubyear
\newcount\@copyear

\bgroup
  \global\@pubyear=\number\year
  \def\tmp@a{\expandafter\@gobbletwo\the\year}
  \global\@copyear\number\tmp@a\relax
\egroup

\def\pubyear#1{\global\@pubyear#1\relax\gdef\@crtyr{#1}%
  \expandafter\global\expandafter\@copyear\@gobbletwo#1\relax
  \ignorespaces}

\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% startlocaldefs, endlocaldefs

\def\startlocaldefs{\makeatletter}
\def\endlocaldefs{\makeatother}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% thanksref, thanksmark, thankslabel, thankstext
\def\myref#1{\expandafter\my@setref\csname r@#1\endcsname\@firstoftwo{#1}}
\let\my@setref\@setref

\long\def\my@footnotetext#1{\insert\footins{%
    \reset@font\footnotesize
    \interlinepenalty\interfootnotelinepenalty
    \splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
    \hsize\columnwidth \@parboxrestore
    \color@begingroup
      \@makefntext{%
        \rule\z@\footnotesep\ignorespaces #1\@finalstrut\strutbox}%
    \color@endgroup}}%


\def\mylabel#1{\@bsphack
  \protected@write\@auxout{}%
         {\string\mynewlabel{#1}{{\@currentlabel}{\thepage}}}%
  \@esphack}

\let\mynewlabel\newlabel

\newcounter{thanks}
\def\thethanks{\arabic{thanks}}

\def\thanksmark@fmt#1{\hbox{$^{#1}$}}
\def\thanksref@sep{,}

% hooks for the hyperref:
\def\thankref@hyperlink#1{\saferef{#1thanks}}
\def\thanks@hypertarget#1{}
\def\thanksmark@link#1#2{#2}


% Isvedame zymes
\DeclareRobustCommand\thanksref{\@ifnextchar[{\@tempswatrue\@thanksref}{\@tempswafalse\@thanksref[]}}

\def\@thanksref[#1]#2{%
  \if@tempswa% []
    \thanksmark@fmt{#1}%
  \else%
    \let\@tempa\@empty%
    \thanksmark@fmt{\@for\@tempb:=#2\do{%
       \@tempa\let\@tempa\thanksref@sep%
       \edef\@tempb{\expandafter\@firstofone\@tempb\@empty}%
       \myref{\@tempb thanks}}}%
   \fi}

% Suformuojame ir isvedame zyme
\def\thanksmark{\@ifnextchar[{\@tempswatrue\@thanksmark}{\@tempswafalse\@thanksmark[]}}

\def\@thanksmark[#1]#2{%
   \if@tempswa% []
     \protected@edef\@currentlabel{#1}%
   \else% 
     \refstepcounter{thanks}%
   \fi%
   \mylabel{#2thanks}%
   \thanksmark@fmt{\expandafter\myref{#2thanks}}}

% Suformuojame tik zyme
\def\thankslabel{\@ifnextchar[{\@tempswatrue\@thankslabel}{\@tempswafalse\@thankslabel[]}}

\def\@thankslabel[#1]#2{%
   \if@tempswa% []
     \protected@edef\@currentlabel{#1}%
   \else% 
     \refstepcounter{thanks}%
   \fi%
   \mylabel{#2thanks}}%

% Suformuojame zyme ir idedame teksta i \@thanks:
\def\thankstext{\@ifnextchar[{\@tempswatrue\@thankstext}{\@tempswafalse\@thankstext[]}}

\def\@thankstext[#1]#2#3{%
  \@thankslabel[#1]{#2}%
  \protected@xdef\@thanks{\@thanks\protect\thanks@thefnmark{#2thanks}\protect\my@footnotetext{#3}}}%

\def\thanks@thefnmark#1{\begingroup\unrestored@protected@xdef\@thefnmark{\myref{#1}}\endgroup}%


% ST makrosas savo numeracijos sistemos sukurimui
\def\setvaluelist#1#2{\@tempcnta=0\relax
  \@for\@curr@val:=#2\do{%
     \advance\@tempcnta by1\relax
     \expandafter\protected@xdef\csname #1@item@\the\@tempcnta\endcsname{\@curr@val}%
     }%
     \expandafter\protected@xdef\csname #1@item@0\endcsname{\the\@tempcnta}%
}
\xdef\getitemvalue#1#2{\noexpand\csname #1@item@#2\endcsname}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% url@fmt

% {url}{text}
\def\ims@href#1#2{#2}

% {prefix}{font}{text}{url}

\def\url@fmt#1#2#3#4{%
   \edef\@tempa{#3}%
   \ifx\@tempa\@empty%
   \else%
     #1{#2\ims@href{#4}{#3}}%
   \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \ead, \printead

\RequirePackage{keyval}

\def\email@text{e-mail: }
\def\url@text{url: }
\def\fullurl@text{url: }
\def\ead@sep{;~}
\let\ead@size\relax
\def\printead@fmt#1{#1}

% naudojame keyval paketa
\define@key{ead}{email}[true]{\def\ead@type{email}}
\define@key{ead}{url}[true]{\def\@tempa{fullurl}\ifx\ead@type\@tempa\else\def\ead@type{url}\fi}
\define@key{ead}{label}{\def\ead@label{#1}}

\define@key{ead}{text}{%
  \bgroup%
    \def\\{\string\break}
    \def\break{\string\break}%
    \protected@edef\@currentlabel{#1}%
    \label{\ead@label @\ead@type text}%
  \egroup}

\define@key{ead}{nopdflink}[true]{%
   \protected@edef\@currentlabel{nolink}%
   \label{\ead@label @nopdflink}}


\DeclareRobustCommand\ead[2][label= ,email]{{%
  \def\ead@type{email}% default
  \checkead@prefix#2://\end%
  \def\texttildelow{\noexpand\texttildelow}%
  \setkeys{ead}{#1}%
  \protected@edef\@currentlabel{#2}%
  \label{\ead@label @\ead@type}}}

\def\checkead@prefix#1://#2\end{\ifx.#2.\else\def\ead@type{fullurl}\fi}

\newif\ifnot@ead@star
\newif\if@printead@opt

\DeclareRobustCommand{\printead}{\@ifstar{\not@ead@starfalse\@printead}{\not@ead@startrue\@printead}}

\def\@printead{\@ifnextchar[{\@printead@opttrue\@@printead}{\@printead@optfalse\@@printead[]}}

\def\@@printead[#1]#2{{%
   \if@printead@opt%[]
      \def\safe@href@text{#1}%
      \not@ead@starfalse%
   \fi%      
   \let\prev@ead@text\relax%
   \let\@ead@sep\relax%
   \let\ead@text\relax%
   \let\ead@prefix\relax%
   \def\ead@type{}%
   \@tempcnta=0%
   \let\sv@safe@href\safe@href%
   \printead@fmt{\@for\ead@ref:=#2\do{%
       \advance\@tempcnta by1%
       \let\safe@href\sv@safe@href%
       \@ead@sep\let\@ead@sep\ead@sep%
       \@ifundefined{r@\ead@ref @nopdflink}{}{\def\safe@href##1##2{##2}}%
       \@ifundefined{r@\ead@ref @email}{}{\let\ead@text\email@text\def\ead@type{email}\def\ead@prefix{mailto:}}%
       \@ifundefined{r@\ead@ref @url}{}{\let\ead@text\url@text\def\ead@type{url}\def\ead@prefix{http://}}%
       \@ifundefined{r@\ead@ref @fullurl}{}{\let\ead@text\fullurl@text\def\ead@type{fullurl}\def\ead@prefix{}}%
       \ifx\prev@ead@text\ead@text\let\ead@text\relax\fi%
       \if@printead@opt\ifnum\@tempcnta>1\@latex@error{Command \@backslashchar printead[]{e1} could have only one parameter "e1"!}\@eha\fi%
       \else\@ifundefined{r@\ead@ref @\ead@type text}{\def\safe@href@text{\saferef{\ead@ref @\ead@type}}}{\def\safe@href@text{\saferef{\ead@ref @\ead@type text}}}\fi%
       \ifnot@ead@star\ead@text\fi{\ead@size\def\null{}\safe@href{\ead@prefix\saferef{\ead@ref @\ead@type}}{\safe@href@text}}%
       \@ifundefined{ead@text}{}{\let\prev@ead@text\ead@text}}}%                        
}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normaltext, nohyphen, today

% normalus tekstas (justify)
\def\normaltext{\let\\=\@normalcr%
  \leftskip\z@ \@rightskip\z@ \rightskip\@rightskip%
  \parfillskip\@flushglue}

% skiemenavimo isjungimas
\def\nohyphen{\pretolerance=\@M \tolerance=\@M \hyphenpenalty=\@M \exhyphenpenalty=\@M}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newpseudoenvironment
% \begin ... \end mechanizmo perrabibrezimas, kad 
% galima butu isvengti papildomu grupiu (ST)

\let\org@begin\begin
\let\org@end\end
\def\begin#1{%
  \@ifundefined{pseudo@#1}%
    {\org@begin{#1}}{\csname pseudo@#1\endcsname[0]\relax}%
  }
\def\end#1{%
  \@ifundefined{pseudo@#1}%
    {\org@end{#1}}{\csname pseudo@#1\endcsname[1]\relax}%
  }
\def\newpseudoenvironment#1#2#3{%
  \expandafter\gdef\csname pseudo@#1\endcsname[##1]{%
     \relax\ifcase##1\relax\def\@@next@@{#2}\or\def\@@next@@{#3}\else\let\@@next@@\relax\fi\@@next@@}%
  }


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dbl@error

\def\dbl@error#1{\@latex@error{#1}\@ehc}

%%%% end   of vtexfm.sty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pakeitimai vtexfmpdf.sty:
\def\saferef#1{\expandafter\safe@setref\csname r@#1\endcsname\@firstoftwo{#1}}
\let\safe@setref\@setref

% {url}{text}
\def\safe@href#1#2{#2}

\def\email@text{E-mail: }
\def\url@text{URL: }
\def\fullurl@text{URL: }


% SECTIONS
\setcounter{secnumdepth}{5}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]

\renewcommand\thesection      {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}

\setcounter{tocdepth}{3}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION commands:
% from latex.ltx:
% Two improvements:
% 1. if section command is defined as "inline" the '.' will be inserted after heading;
% 2. section* will write to toc and will appear in pdf bookmarks

% dirty trick...
\def\@startsection#1#2#3#4#5#6{%
  \if@noskipsec \leavevmode \fi
  \par
  \@tempskipa #4\relax
  \@afterindenttrue
  \ifdim \@tempskipa <\z@
    \@tempskipa -\@tempskipa \@afterindentfalse
  \fi
  \if@nobreak
    \everypar{}%
  \else
    \addpenalty\@secpenalty\addvspace\@tempskipa
  \fi
  \@ifstar
    {\def\ssection@level{#2}\@ssect{#3}{#4}{#5}{#6}}%
    {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}

% Trick for the hyperref:
\def\setaftersec@dot#1{\if@autosecdot\setbox0=\hbox{#1}\ifdim\wd0>0\p@\def\aftersec@dot{.}\else\def\aftersec@dot{}\fi\fi}
\let\aftersec@dot\relax

% section - will add hook for the dot after section heading
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
    \protected@edef\@svsec{\@seccntformat{#1}\relax}%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #6{%
        \@hangfrom{\hskip #3\relax\@svsec}%
          \interlinepenalty \@M #8\@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7}%
  \else
    \setaftersec@dot{#8}%
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec #8\aftersec@dot}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}

% section* - will add hook for the dot after section heading and \contentsline 
\def\@ssect#1#2#3#4#5{%
  \@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M #5\@@par}%
    \endgroup
  \else
    \setaftersec@dot{#5}%
    \def\@svsechd{#4{\hskip #1\relax #5\aftersec@dot}}%
  \fi
  \ifnum\ssection@level=1\phantomsection\addcontentsline{toc}{section}{#5}\fi%
  \@xsect{#3}}

% Block adding to contents for the next command only:
\def\nocontentsline{%
  \let\@@addcontentsline\addcontentsline%
  \def\addcontentsline##1##2##3{\let\addcontentsline\@@addcontentsline}}

% \phantomsection is defined in hyperref
\let\phantomsection\relax



%%%% start of vtexlst.sty [2002/02/09 VTEX list type environments (VS)]

%%%%%%%%%%%%%% Hooks
% general list indent
\def\list@parindent{1pc}

% enumerate
\def\enumerate@indent{\z@}
\def\enumerate@parindent{\list@parindent}
\def\enumerate@labelsep{.5em}
\def\enumerate@firsttopsep{\smallskipamount}
\def\enumerate@topsep{\z@}
\def\enumerate@rightskip{\z@}% {\z@ plus 1fill}
\let\enumerate@hook\relax

% longlist
\def\longlist@indent{\list@parindent}
\def\longlist@parindent{\list@parindent}
\def\longlist@labelsep{.5em}
\def\longlist@topsep{\smallskipamount}
\def\longlist@rightskip{\z@}% {\z@ plus 1fill}
\let\longlist@hook\relax

% itemize
\def\itemize@indent{\z@}
\def\itemize@parindent{\list@parindent}
\def\itemize@labelsep{.5em}
\def\itemize@firsttopsep{\smallskipamount}
\def\itemize@topsep{\z@}
\def\itemize@rightskip{\z@}% {\z@ plus 1fill}
\let\itemize@hook\relax

% description
\def\description@indent{\list@parindent}
\let\description@topsep\smallskipamount
\let\description@hook\relax

% quotation
\let\quotation@size\footnotesize
\def\quotation@itemindent{\list@parindent}
\def\quotation@parindent{\list@parindent}
\def\quotation@leftmargin{\list@parindent}
\let\quotation@rightmargin\z@
\let\quotation@topsep\smallskipamount


% quote
\let\quote@size\footnotesize
\def\quote@indent{\z@}
\def\quote@leftmargin{2pc}
\def\quote@rightmargin{2pc}
\let\quote@topsep\smallskipamount

% verse
\def\verse@margin{2pc}
\let\verse@topsep\smallskipamount


%%%%%%%%%%%%%% Parameters
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty


%%%%%%%%%%%%%% General list parameters (not important}

\parsep\z@
\topsep\z@
\partopsep\z@
\itemsep\z@
\labelsep.5em

\def\@listi  {\leftmargin\leftmargini}

\def\@listii {\leftmargin\leftmarginii%
             \labelwidth\leftmarginii%
             \advance\labelwidth-\labelsep}%

\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep}

\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}

\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}

\setlength\leftmargini   {2em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv   {1em}
\setlength\leftmarginvi  {1em}

\@listi

\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}

%%%%%%%%%%%%%% Enumerate

\def\theenumi  {\@arabic\c@enumi}
\def\theenumii {\@alph\c@enumii}
\def\theenumiii{\@roman\c@enumiii}
\def\theenumiv {\@Alph\c@enumiv}

\def\labelenumi  {\theenumi.}
\def\labelenumii {\theenumii)}
\def\labelenumiii{\theenumiii.}
\def\labelenumiv {\theenumiv.}

\def\enumerate{\@ifnextchar[{\@enumerate}{\@enumerate[]}}

\def\@enumerate[#1]{%
    \ifnum \@enumdepth >\thr@@ \@toodeep\else
    \advance\@enumdepth \@ne 
    \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
    \list {\csname label\@enumctr\endcsname}%
       {\usecounter{\@enumctr}
        \def\makelabel##1{{\hfill\normalfont ##1}}
%
        \ifx.#1.%
           \setcounter{\@enumctr}{3}% 
           \expandafter\settowidth\expandafter\labelwidth\expandafter{\csname label\@enumctr\endcsname}%
           \setcounter{\@enumctr}{0}%
        \else%
           \settowidth{\labelwidth}{#1}
        \fi%
        \labelsep=\enumerate@labelsep%
        \itemindent=\labelwidth \advance\itemindent\labelsep%
        \ifnum \@enumdepth=1\advance\itemindent by\enumerate@indent\fi%
        \leftmargin=\the\itemindent%
        \itemindent\z@%
        \partopsep\z@  \parsep\z@ \itemsep\z@ \rightmargin\z@
        \ifnum \@enumdepth=1 \topsep\enumerate@firsttopsep \else \topsep\enumerate@topsep\fi%
        \listparindent\enumerate@parindent%        
        \@rightskip\enumerate@rightskip%
        \enumerate@hook% 
        }\fi}

\let\endenumerate =\endlist

%%%%%%%%%%%%%% Longlist

\newcounter{longlist}
\def\thelonglist   {\@roman\c@longlist}
\def\labellonglist {\thelonglist)}

\def\longlist{\@ifnextchar[{\@longlist}{\@longlist[]}}

\def\@longlist[#1]{%
    \list {\labellonglist}%
       {\usecounter{longlist}
        \def\makelabel##1{{\hfill\normalfont ##1}}
%
        \ifx.#1.%
           \setcounter{longlist}{3}% 
           \settowidth\labelwidth{\labellonglist}%
           \setcounter{longlist}{0}%
        \else%
           \settowidth{\labelwidth}{#1}
        \fi%
        \labelsep=\longlist@labelsep%
        \itemindent=\labelwidth \advance\itemindent\labelsep%
        \advance\itemindent by\longlist@indent%
        \rightmargin\z@ \leftmargin\z@%
        \partopsep\z@  \parsep\z@ \itemsep\z@ %
        \topsep\longlist@topsep%
        \listparindent\longlist@parindent%        
        \@rightskip\longlist@rightskip%
        \longlist@hook% 
        }}

\let\endlonglist =\endlist


%%%%%%%%%%%%%%  setenumtype

\def\setenumtype#1{%
\if#1a \let\enumerate\longlist\fi%
\if#1c \def\enumerate@indent{\z@}\fi}%

%%%%%%%%%%%%%%  itemize

\def\labelitemi  {\textbullet}
\def\labelitemii {\textendash}
\def\labelitemiii{\textasteriskcentered} 
\def\labelitemiv {{\footnotesize +}}

\def\itemize{\@ifnextchar[{\@itemize}{\@itemize[]}}

\def\@itemize[#1]{%
    \ifnum \@itemdepth >\thr@@ \@toodeep\else
    \advance\@itemdepth \@ne 
    \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
    \list {\csname\@itemitem\endcsname}%
       {\def\makelabel##1{{\hfill\normalfont ##1}}        
%
        \ifx.#1.%
           \settowidth\labelwidth{\csname\@itemitem\endcsname}%
        \else%
           \settowidth{\labelwidth}{#1}
        \fi%
        \labelsep=\itemize@labelsep%
        \itemindent=\labelwidth \advance\itemindent\labelsep%
        \ifnum \@itemdepth=1\advance\itemindent by\itemize@indent\fi%
        \leftmargin=\the\itemindent%
        \itemindent\z@%
        \partopsep\z@  \parsep\z@ \itemsep\z@ \rightmargin\z@
        \ifnum \@itemdepth=1 \topsep\itemize@firsttopsep \else \topsep\itemize@topsep\fi%
        \listparindent\itemize@parindent%        
        \@rightskip\itemize@rightskip%
        \itemize@hook% 
        }\fi}

\let\enditemize =\endlist

%%%%%%%%%%%%%%  description

\def\description{%
         \list{}{\labelwidth\z@ \labelsep\z@  \rightmargin\z@
         \leftmargin\z@ \advance\leftmargin by\description@indent
         \itemindent-\description@indent
         \partopsep\z@ \topsep\description@topsep \parsep\z@
         \itemsep\z@ \listparindent\z@
         \description@hook%
         \let\makelabel\descriptionlabel}}

\def\enddescription{\endlist}

\def\descriptionlabel#1{\if.#1.\else\bfseries #1\hskip.5em\fi}

%%%%%%%%%%%%%%  quotation

\def\quotation{%
        \list{}{\quotation@size%
        \listparindent\quotation@parindent%
        \itemindent   \quotation@itemindent%
        \rightmargin\quotation@rightmargin   \leftmargin\quotation@leftmargin%
        \partopsep\z@ \topsep\quotation@topsep \parsep\z@%
                        }%
        \item[\Q@strut]\relax}

\def\endquotation{\endlist}

\def\Q@strut{\leavevmode\hbox{\vrule height9pt depth1pt width0pt}}

%%%%%%%%%%%%%%  quote

\def\quote{%
        \list{}{\quote@size%
        \listparindent\quote@indent%
        \itemindent \listparindent%
        \rightmargin\quote@rightmargin   \leftmargin\quote@leftmargin%
        \partopsep\z@ \topsep\quote@topsep \parsep\z@%
                       }%
        \item\relax}

\def\endquote{\endlist}

%%%%%%%%%%%%%%  verse


\def\verse{%
   \let\\\@centercr
   \list{}{\itemindent\z@ \listparindent\z@ \rightmargin\verse@margin \leftmargin\rightmargin
	   \partopsep\z@ \parsep\z@ \topsep\verse@topsep \itemsep\z@}%
   \item\relax}

\def\endverse{\endlist}

%%%% end   of vtexlst.sty
\renewcommand\labellonglist {(\thelonglist)}

%%%% start of vtexmth.sty [2003/05/08 VTEX math settings (VS)]
%%%%%%%%%%%%%% Parameters
\setlength\arraycolsep   {1.2\p@}
\setlength\tabcolsep     {6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep {2\p@}
\setlength\tabbingsep    {6\p@}
\setlength\jot           {6\p@}


%%%%%%%%%%%%%% Math spacing
\setlength\thinmuskip  {2.5mu}
\setlength\medmuskip   {3.5mu plus 1.5mu minus 2mu}
\setlength\thickmuskip {3.75mu plus 1.5mu}
\setlength\mathsurround{0.1\p@}

%%%%%%%%%%%%%% Equation numbering

\ifcase\eqn@numbering
   \renewcommand\theequation{\@arabic\c@equation}
\or %seceqn  
  \@addtoreset {equation}{section}
  \renewcommand\theequation{\thesection.\@arabic\c@equation}
\or %chapeqn
  \@addtoreset {equation}{chapter}
  \renewcommand\theequation{\thechapter.\@arabic\c@equation}
\fi

%%%%%%%%%%%%%% tfrac, dfrac, dvt, dvtx

\def\tfrac#1#2{{\textstyle\frac{#1}{#2}}}
\def\dfrac#1#2{{\displaystyle\frac{#1}{#2}}}

\def\dvt{\colon\ }

\def\dvtx{%
  \mathchoice{\nobreak\,\colon\relax}%
             {\nobreak\,\colon\relax}%
             {\nobreak\,\colon\;\relax}%
             {\nobreak\,\colon\;\relax}%
  }


%%%%%%%%%%%%%% tsub, textsubscript, tsup, textsuperscript

\let\textsuperscript\relax
\DeclareRobustCommand*\textsuperscript[1]{%
  \@textsuperscript{\selectfont#1}}
\def\@textsuperscript#1{%
  {\m@th\ensuremath{^{\mathchoice{\mbox{\fontsize\f@size\z@#1}}%
                                 {\mbox{\fontsize\tf@size\z@#1}}%
                                 {\mbox{\fontsize\sf@size\z@#1}}%
                                 {\mbox{\fontsize\ssf@size\z@#1}}%
                     }%
                   }%
  }%
}

\DeclareRobustCommand*\textsubscript[1]{%
  \@textsubscript{\selectfont#1}}
\def\@textsubscript#1{%
  {\m@th\ensuremath{_{\mathchoice{\mbox{\fontsize\f@size\z@#1}}%
                                 {\mbox{\fontsize\tf@size\z@#1}}%
                                 {\mbox{\fontsize\sf@size\z@#1}}%
                                 {\mbox{\fontsize\ssf@size\z@#1}}%
                     }%
                   }%
  }%
}


\let\tsub\textsubscript
\let\tsup\textsuperscript

%%%%%%%%%%%%%% big, bigg

\newcommand{\m@kebig}[3]{\hbox{#1$\left#2\vbox to #3ex{}\right.\n@space$}}
\newcommand{\m@kelarge}[2]{\mathchoice{\m@kebig{}{#1}{#2}}{\m@kebig{}{#1}{#2}}
    {\m@kebig\scriptsize{#1}{#2}}{\m@kebig\tiny{#1}{#2}}}
\renewcommand{\big}[1]{{\m@kelarge{#1}{1.79999}}}
\renewcommand{\Big}[1]{{\m@kelarge{#1}{2.43529}}}
\renewcommand{\bigg}[1]{{\m@kelarge{#1}{3.07059}}}
\renewcommand{\Bigg}[1]{{\m@kelarge{#1}{3.70588}}}

%%%%%%%%%%%%%% divid, zs

\def\divid#1#2{%
  \mathchoice{{#1}\over{#2}}% Display mode
    {\vbox{\hbox{$\scriptstyle #1$}\vskip1pt}\over
       \vbox{\vskip1pt\hbox{$\scriptstyle #2$}\vskip1pt}}% Text mode
    {\vbox{\hbox{$\scriptstyle #1$}\vskip1pt}\over
       \vbox{\vskip2pt\hbox{$\scriptstyle #2$}\vskip1pt}}% Script mode
    {\vbox{\hbox{$\scriptstyle #1$}\vskip1pt}\over
       \vbox{\vskip.7pt\hbox{$\scriptstyle #2$}\vskip3pt}}% Sriptsript mode
 }
% ZS apibrezimas
\def\zs#1{\mathchoice{#1}{#1}{\lower.4ex\hbox{$\scriptstyle#1$}}
   {\lower0.4ex\hbox{$\scriptscriptstyle#1$}}}
\def\zss#1{\mathchoice{#1}{#1}{\lower.8ex\hbox{$\scriptstyle#1$}}
   {\lower0.8ex\hbox{$\scriptscriptstyle#1$}}}

%%%%%%%%%%%%%% operatorname
\def\operatorname#1{\mathop{\mathrm{#1}}\nolimits}

%%%%%%%%%%%%%% bolds

\DeclareRobustCommand{\bolds}[1]{%
  \begingroup
  \let\@nomath\@gobble \mathversion{bold}%
  \math@atom{#1}{%
  \mathchoice%
    {\hbox{$\m@th\displaystyle#1$}}%
    {\hbox{$\m@th\textstyle#1$}}%
    {\hbox{$\m@th\scriptstyle#1$}}%
    {\hbox{$\m@th\scriptscriptstyle#1$}}}%
  \endgroup}

\def\math@atom#1#2{%
   \binrel@{#1}\binrel@@{#2}}

\def\binrel@#1{\begingroup
  \setbox\z@\hbox{\thinmuskip0mu
    \medmuskip\m@ne mu\thickmuskip\@ne mu
    \setbox\tw@\hbox{$#1\m@th$}\kern-\wd\tw@
    ${}#1{}\m@th$}%
  \edef\@tempa{\endgroup\let\noexpand\binrel@@
    \ifdim\wd\z@<\z@ \mathbin
    \else\ifdim\wd\z@>\z@ \mathrel
    \else \relax\fi\fi}%
  \@tempa
  }
\let\binrel@@\relax
\let\boldsymbol\bolds

%%%%%%%%%%%%%% tabular

\def\ccline#1{\do@ccline(#1,1-0)}
\def\do@ccline(#1-#2,#3){%
  \noalign{\ifnum#1>#2\gdef\temp@next{}\else\gdef\temp@next{\@cline#1-#2\@nil\do@ccline(#3,1-0)}\fi}%
  \temp@next
  }


%%%%%%%%%%%%%% aligned

\newenvironment{aligned}{%
  \relax\ifmmode\else{\@latex@error{aligned aplinka naudojama tik matematinej modoje!}\@ehc}\fi
  \intertext@\null\@ifnextchar[{\aligned@}{\aligned@[c]}%
  }{\crcr\egroup\egroup}

\def\aligned@[#1]{%
  \if #1t\vtop \else \if#1b\vbox \else \vcenter \fi\fi
  \bgroup
   \let\\\c@eqncr
   \ialign\bgroup\hfil\strut$\m@th\displaystyle{##{}}$&&%
    $\m@th\displaystyle{{}##{}}$\hfil\crcr}

\def\c@eqncr{%
  {\ifnum0=`}\fi
   \@ifstar{\global\@eqpen\@M\c@yeqncr
     }{\global\@eqpen\interdisplaylinepenalty \c@yeqncr}%
  }

\def\c@yeqncr{\@testopt\c@xeqncr\z@skip}
\def\c@xeqncr[#1]{%
  \ifnum0=`{\fi}%
  \cr\noalign{\penalty\@eqpen\vskip\jot\vskip #1\relax}%
  }

\def\intertext{\ClassError{fied}{Invalid use of \string\intertext}\@eha}
\def\intertext@{\def\intertext##1{\noalign{%
  \penalty\postdisplaypenalty\vskip\belowdisplayskip
  \vbox{\normalbaselines\noindent##1}%
  \penalty\predisplaypenalty\vskip\abovedisplayskip\vskip\jot}}}



%%%% end   of vtexmth.sty


% problema: neteisingoje vietoje numeris, kai formule list'e:
\renewcommand\@eqnnum{\hb@xt@.01\p@{}%
                      \rlap{\normalfont\normalcolor
                       \hskip-\displayindent\hskip -\displaywidth(\theequation)}}


% triukas ikisti teksta i numeri:
\def\eqntext#1{\global\let\curr@eqnnum\@eqnnum\gdef\@eqnnum##1##2{$#1$\global\let\@eqnnum\curr@eqnnum}}

%%%% start of vtexthm.sty [2001/11/24 VTEX theorem type environments (VS)]

%%%%%%%%%%%%%% HOOKS
  \newskip\theorempreskipamount
  \newskip\theorempostskipamount
  \global\setlength\theorempreskipamount{12\p@ \@plus 2\p@ \@minus 1\p@}
  \global\setlength\theorempostskipamount{12\p@ \@plus 2\p@ \@minus 1\p@}
  \let\thm@indent\z@%
  \def\@thmcounter#1{\noexpand\arabic{#1}}

  \gdef\thmcounter@sep{.} % taskas po teoremos skaitliuko
  \gdef\xxthmcounter@sep{.} % taskas po teoremos skaitliuko xxtheorem atveju

  \let\thm@skip\enskip   % tarpas nuo teoremos iki teksto

% Proof:
  \def\proofname{Proof}
  \def\proof@sep{.}
  \let\proof@shape\scshape
  \let\proof@skip\enskip
  \let\proof@indent\noindent

% NEWTHEOREM and PROCLAIM environments 
%
\gdef\theoremstyle#1{%
   \@ifundefined{th@#1}{\@warning{Unknown theoremstyle `#1'. Using `plaindot'}%
          \theorem@style{plaindot}}%
      {\theorem@style{#1}}%
      \begingroup
        \csname th@\the\theorem@style \endcsname
      \endgroup}

\newtoks\theorem@style
\global\theorem@style{plaindot}

\gdef\@xnthm#1#2[#3]{\expandafter\@ifdefinable\csname #1\endcsname
   {%
    \@definecounter{#1}\@addtoreset{#1}{#3}%
    \expandafter\xdef\csname the#1\endcsname
      {\expandafter \noexpand \csname the#3\endcsname
       \thmcounter@sep \@thmcounter{#1}}%
    \def\@tempa{\global\@namedef{#1}}%
    \expandafter \@tempa \expandafter{%
      \csname th@\the \theorem@style \endcsname
     \@thm{#1}{#2}}%
    \global \expandafter \let \csname end#1\endcsname \@endtheorem
    \def\@tempb{\global\@namedef{#1*}}%
    \expandafter \@tempb \expandafter{%
      \csname th@\the \theorem@style \endcsname
     \@thmstar{#2}}%
    \global \expandafter \let \csname end#1*\endcsname \@endtheorem
   }}

\gdef\@ynthm#1#2{\expandafter\@ifdefinable\csname #1\endcsname
   {\@definecounter{#1}%
    \expandafter\xdef\csname the#1\endcsname{\@thmcounter{#1}}%
    \def\@tempa{\global\@namedef{#1}}\expandafter \@tempa
     \expandafter{\csname th@\the \theorem@style \endcsname
    \@thm{#1}{#2}}%
    \global \expandafter \let \csname end#1\endcsname \@endtheorem
    \def\@tempb{\global\@namedef{#1*}}\expandafter \@tempb
     \expandafter{\csname th@\the \theorem@style \endcsname
    \@thmstar{#2}}%
    \global \expandafter \let \csname end#1*\endcsname \@endtheorem}}

\gdef\@othm#1[#2]#3{\expandafter\@ifdefinable\csname #1\endcsname
  {\expandafter \xdef \csname the#1\endcsname
     {\expandafter \noexpand \csname the#2\endcsname}%
    \def\@tempa{\global\@namedef{#1}}\expandafter \@tempa
     \expandafter{\csname th@\the \theorem@style \endcsname
    \@thm{#2}{#3}}%
    \global \expandafter \let \csname end#1\endcsname \@endtheorem
    \def\@tempb{\global\@namedef{#1*}}\expandafter \@tempb
     \expandafter{\csname th@\the \theorem@style \endcsname
    \@thmstar{#3}}%
    \global \expandafter \let \csname end#1*\endcsname \@endtheorem}}

\gdef\@thm#1#2{\refstepcounter{#1}%
   \trivlist
   \itemindent\thm@indent%
   \@topsep \theorempreskipamount%
   \@topsepadd \theorempostskipamount%
   \@ifnextchar [%
   {\@ythm{#1}{#2}}%
   {\@begintheorem{#2}{\csname the#1\endcsname}\ignorespaces}}

\gdef\@thmstar#1{%
   \trivlist
   \itemindent\thm@indent%
   \@topsep \theorempreskipamount%
   \@topsepadd \theorempostskipamount%
   \@ifnextchar [%
   {\@ythmstar{#1}}%
   {\@@begintheorem{#1}\ignorespaces}}

\def\@ythmstar#1[#2]{%
  \@@opargbegintheorem{#1}{#2}\ignorespaces}

\global\let\@xthm\relax

\if@xx@theorem
  \def\fmt@thm@name#1#2{#2\xxthmcounter@sep\ #1}
\else
  \def\fmt@thm@name#1#2{#1\ #2}
\fi

\def\theorem@@shape#1#2#3{\expandafter\def\csname th@#1@#2@shape\endcsname{#3}}
\let\theoremshape\theorem@@shape

\theorem@@shape{plaindot}{head}{\normalfont\scshape}
\theorem@@shape{plaindot}{optional}{\mdseries\upshape}
\theorem@@shape{plaindot}{body}{\normalfont\itshape}

\theorem@@shape{proclaim}{head}{\normalfont\scshape}
\theorem@@shape{proclaim}{optional}{\normalfont}
\theorem@@shape{proclaim}{body}{\normalfont\upshape}


\gdef\th@plaindot{%
  \let\th@main@head@shape\th@plaindot@head@shape
  \let\th@main@optional@shape\th@plaindot@optional@shape
  \let\th@main@body@shape\th@plaindot@body@shape
  \th@main@theorem
}

\gdef\th@proclaim{%
  \let\th@main@head@shape\th@proclaim@head@shape
  \let\th@main@optional@shape\th@proclaim@optional@shape
  \let\th@main@body@shape\th@proclaim@body@shape
  \th@main@theorem
}

\gdef\th@main@theorem{\th@main@body@shape
  \def\@begintheorem##1##2{%
        \item\relax{\normalfont\th@main@head@shape\fmt@thm@name{##1}{##2}\thmcounter@sep}\thm@skip}%
  \def\@@begintheorem##1{%
        \item\relax{\normalfont\th@main@head@shape ##1\thmcounter@sep}\thm@skip}%
\def\@opargbegintheorem##1##2##3{%
   \item\relax{\normalfont\th@main@head@shape\fmt@thm@name{##1}{##2}}\ {\th@main@optional@shape ##3\thmcounter@sep}\thm@skip}%
\def\@@opargbegintheorem##1##2{%
   \item\relax{\normalfont\th@main@head@shape ##1}\ {\th@main@optional@shape ##2\thmcounter@sep}\thm@skip}}

\let\old@newtheorem\newtheorem
\def\newtheorem{\theoremstyle{plaindot}\old@newtheorem}
\def\newproclaim{\theoremstyle{proclaim}\old@newtheorem}

\ifcase\thm@numbering
   \newtheorem{thm}{Theorem}
   \def\thethm{\@arabic\c@thm}
\or % secthm
  \newtheorem{thm}{Theorem}[section]
  \@addtoreset{thm}{section}
  \def\thethm{\thesection.\@arabic\c@thm}
\or % chapthm
  \newtheorem{thm}{Theorem}[chapter]
  \@addtoreset{thm}{chapter}
  \def\thethm{\thechapter.\@arabic\c@thm}
\fi

% PROOF ENVIRONMENT

\def\theqed{\ensuremath{\square}}

\gdef\qed{\relax\ifmmode~\hfill\theqed\else\unskip\nobreak~\hfill\theqed\fi}
\let\@qed\qed

\gdef\qedbreak{\par\nobreak\parfillskip0pt\hfill}

\gdef\upqed{\vskip-\belowdisplayskip\vskip-\baselineskip\kern0pt}

\def\noqed{\let\sv@qed\@qed\def\@qed{\global\let\@qed\sv@qed}}

\newenvironment{pf}%
  {\par\addvspace{\theorempreskipamount}\proof@indent
   {\proof@shape\proofname\proof@sep}\proof@skip\ignorespaces}%
   {\@qed\par\addvspace{\theorempreskipamount}}

\@namedef{pf*}#1{\par\begingroup\def\proofname{#1}\pf\endgroup\ignorespaces}
\@namedef{endpf*}{\@qed\par\addvspace{\theorempostskipamount}}

%%%% end   of vtexthm.sty

  \setattribute{proof}{indent}{}
  \setattribute{thm}{indent}{\parindent}
  \setattribute{thm}{skip}{\hskip1em}
  \setattribute{proof}{skip}{\hskip1em}
  \gdef\qed{\relax\ifmmode~\hfill\square\else\unskip\hskip10pt$\square$\fi}
  \let\@qed\qed
  \def\mathend@qed{\hbox to\textwidth{\hfill$\square$}}
  \gdef\upqed{\vskip-\belowdisplayskip\vskip-\baselineskip\kern0pt\let\@qed\mathend@qed}
  \gdef\rightqed{\parfillskip0pt\par\nobreak\let\@qed\mathend@qed}

  \global\setlength\theorempreskipamount\medskipamount
  \global\setlength\theorempostskipamount\medskipamount

  \theorem@@shape{plaindot}{body}{\normalfont\itshape}
  \def\@endtheorem{\endtrivlist\@endpefalse}


% BIBLIOGRAPHY
  \def\thebibliography{}
  \def\newblock{\hskip .11em\@plus.33em\@minus.07em}

  \if@nameyear
      \RequirePackage[comma]{natbib}
  \else
    \RequirePackage[square,comma,numbers,openbib,sort&compress]{natbib}
    \newlength{\bibindent}
    \setlength{\bibindent}{20pt}
    \let\cite\citep

  \fi

   \def\bibsection{\specialsection*{\bibname}}
   \let\bibfont\footnotesize
   \def\authors#1{#1}

%%%%%%%%%% change to natbib
% we will make links only from year - will solve problem with 
% author names sticking in the margin

\def\NAT@citex%
  [#1][#2]#3{%
  \NAT@reset@parser
  \NAT@sort@cites{#3}%
  \NAT@reset@citea
  \@cite{\let\NAT@nm\@empty\let\NAT@year\@empty
    \@for\@citeb:=\NAT@cite@list\do
    {\@safe@activestrue
     \edef\@citeb{\expandafter\@firstofone\@citeb\@empty}%
     \@safe@activesfalse
     \@ifundefined{b@\@citeb\@extra@b@citeb}{\@citea%
       {\reset@font\bfseries ?}\NAT@citeundefined
                 \PackageWarning{natbib}%
       {Citation `\@citeb' on page \thepage \space undefined}\def\NAT@date{}}%
     {\let\NAT@last@nm=\NAT@nm\let\NAT@last@yr=\NAT@year
      \NAT@parse{\@citeb}%
      \ifNAT@longnames\@ifundefined{bv@\@citeb\@extra@b@citeb}{%
        \let\NAT@name=\NAT@all@names
        \global\@namedef{bv@\@citeb\@extra@b@citeb}{}}{}%
      \fi
     \ifNAT@full\let\NAT@nm\NAT@all@names\else
       \let\NAT@nm\NAT@name\fi
     \ifNAT@swa\ifcase\NAT@ctype
       \if\relax\NAT@date\relax
         \@citea\NAT@hyper@{\NAT@nmfmt{\NAT@nm}\NAT@date}%
       \else
         \ifx\NAT@last@nm\NAT@nm\NAT@yrsep
            \ifx\NAT@last@yr\NAT@year
              \def\NAT@temp{{?}}%
              \ifx\NAT@temp\NAT@exlab\PackageWarningNoLine{natbib}%
               {Multiple citation on page \thepage: same authors and
               year\MessageBreak without distinguishing extra
               letter,\MessageBreak appears as question mark}\fi
              \NAT@hyper@{\NAT@exlab}%
            \else\unskip\NAT@spacechar
              \NAT@hyper@{\NAT@date}%
            \fi
         \else
           \@citea%
             \NAT@nmfmt{\NAT@nm}%
             \NAT@aysep\NAT@spacechar%
             \NAT@hyper@{\NAT@date}%
         \fi
       \fi
     \or\@citea\NAT@hyper@{\NAT@nmfmt{\NAT@nm}}%
     \or\@citea\NAT@hyper@{\NAT@date}%
     \or\@citea\NAT@hyper@{\NAT@alias}%
     \fi \NAT@def@citea
     \else
       \ifcase\NAT@ctype
        \if\relax\NAT@date\relax
          \@citea\NAT@hyper@{\NAT@nmfmt{\NAT@nm}}%
        \else
         \ifx\NAT@last@nm\NAT@nm\NAT@yrsep
            \ifx\NAT@last@yr\NAT@year
              \def\NAT@temp{{?}}%
              \ifx\NAT@temp\NAT@exlab\PackageWarningNoLine{natbib}%
               {Multiple citation on page \thepage: same authors and
               year\MessageBreak without distinguishing extra
               letter,\MessageBreak appears as question mark}\fi
              \NAT@hyper@{\NAT@exlab}%
            \else
              \unskip\NAT@spacechar
              \NAT@hyper@{\NAT@date}%
            \fi
         \else
           \@citea%
             \NAT@nmfmt{\NAT@nm}%
             \NAT@spacechar\NAT@@open\if*#1*\else#1\NAT@spacechar\fi%
             \NAT@hyper@{\NAT@date}%
         \fi
        \fi
       \or\@citea\NAT@hyper@{\NAT@nmfmt{\NAT@nm}}%
       \or\@citea\NAT@hyper@{\NAT@date}%
       \or\@citea\NAT@hyper@{\NAT@alias}%
       \fi
       \if\relax\NAT@date\relax
         \NAT@def@citea
       \else
         \NAT@def@citea@close
       \fi
     \fi
     }}\ifNAT@swa\else\if*#2*\else\NAT@cmt#2\fi
     \if\relax\NAT@date\relax\else\NAT@@close\fi\fi}{#1}{#2}}



%%%%%%%%%%

%% Load hyperref:

\if@loadhyperref
   \AtEndOfClass{\RequirePackage[dvips]{hyperref}}
%   \AtEndOfClass{\RequirePackage{hypdvips}}
\fi

% arXiv hyperref versija sena, nepalaiko \texttildelow
\DeclareTextSymbol{\texttildelow}{OT1}{126}


% FLOATS: 


%%%% start of vtexflt.sty [2002/01/18 VTEX float type environments (VS)]

\def\figurename{Figure}
\def\tablename {Table}

% Figure
\let\figurefloat@shape\centering
\def\figurecaption@size{\footnotesize}
\def\figurecaption@skip{8\p@}
\def\figurename@skip{\hskip.5em}
\let\figurename@size\relax

% Table
\let\starttable@hook\relax
\let\tablefloat@shape\centering
\let\table@size\footnotesize
\def\tablecaption@shape{\centering}%
\def\tablecaption@size{\footnotesize}
\def\tablecaption@skip{4\p@}
\def\tablename@skip{\hskip.5em}
\let\tablename@size\relax

% Tabnote
\def\firsttabnote@skip{\vskip4\p@}
\def\tabnote@skip{\lineskip2\p@}
\let\tabnote@size\footnotesize
\let\tabnote@indent\leavevmode

%%%%%%%%%%%%%% FLOAT PARAMETERS

\setlength\floatsep       {18\p@ \@plus 6\p@ \@minus 2\p@}
\setlength\textfloatsep   {18\p@ \@plus 6\p@ \@minus 2\p@}
\setlength\intextsep      {18\p@ \@plus 6\p@ \@minus 2\p@}
\setlength\dblfloatsep    {18\p@ \@plus 6\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{24\p@ \@plus 6\p@ \@minus 2\p@}

%  For floats on a separate float page or column:
\setlength\@fptop   {0\p@ \@plus 1fil}
\setlength\@fpsep   {8\p@ \@plus 1000fil}
\setlength\@fpbot   {0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 1000fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}

%  Float counters:
\setcounter{topnumber}{5}
\renewcommand\topfraction{.99}
\setcounter{bottomnumber}{5}
\renewcommand\bottomfraction{.99}
\setcounter{totalnumber}{10}
\renewcommand\textfraction{.05}
\renewcommand\floatpagefraction{.9}
\setcounter{dbltopnumber}{5}
\renewcommand\dbltopfraction{.99}
\renewcommand\dblfloatpagefraction{.8}


%%%%%%%%%%%%%%  FIGURES

\newcount\c@figure
\let\cl@figure\@empty
\let\p@figure\@empty

\ifcase\float@numbering
   \def\thefigure{\@arabic\c@figure}
\or %secfloat
   \@addtoreset{figure}{section}
   \def\thefigure{\thesection.\@arabic\c@figure}
\or %chapfloat
   \@addtoreset{figure}{chapter}
   \def\thefigure{\thechapter.\@arabic\c@figure} 
\fi

\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure.}

\def\figure{\let\@makecaption\@makefigurecaption%
            \let\@floatboxreset\figure@boxreset%
            \@float{figure}}

\def\endfigure{\end@float}

\@namedef{figure*}{\let\@makecaption\@makefigurecaption%
                   \let\@floatboxreset\figure@boxreset%
                   \@dblfloat{figure}}

\@namedef{endfigure*}{\end@dblfloat}

\def\figure@boxreset {%
        \reset@font
        \normalsize
        \@setminipage
        \figurefloat@shape}

\long\def\@makefigurecaption#1#2{%
  \vskip\figurecaption@skip%
  \figurecaption@size%
  \setbox\@tempboxa\hbox{{\figurename@size #1}\figurename@skip #2}%
  \ifdim \wd\@tempboxa >\hsize
    \unhbox\@tempboxa\par
  \else
     \figurecaption@fmt{\box\@tempboxa}\par
  \fi}

\def\figurecaption@fmt#1{\hbox to\hsize{\hfill#1\hfill}}

%%%%%%%%%%%%%%  TABLES

\newcount\c@table
\let\cl@table\@empty
\let\p@table\@empty

\ifcase\float@numbering
    \def\thetable{\@arabic\c@table}
\or %secfloat
  \@addtoreset{table}{section}
  \def\thetable{\thesection.\@arabic\c@table}
\or %chapfloat
  \@addtoreset{table}{chapter}
  \def\thetable{\thechapter.\@arabic\c@table}
\fi

\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable.}

\def\table{\let\@makecaption\@maketablecaption\starttable@hook%
           \tablewidth\hsize%
           \let\@floatboxreset\table@boxreset\@float{table}}

\def\endtable{\end@float}

\@namedef{table*}{\let\@makecaption\@maketablecaption\starttable@hook%
                  \tablewidth\hsize%
                  \let\@floatboxreset\table@boxreset\@dblfloat{table}}

\@namedef{table*}{\end@dblfloat}

\def\table@boxreset{%
        \reset@font%
        \tablefloat@shape%
        \table@size
        \@setnobreak
        \@setminipage
        \setcounter{tabnote}{0}}

\newdimen\tablewidth \tablewidth\textwidth

\long\def\@maketablecaption#1#2{
    \tablecaption@shape%
    \parbox[t]{\tablewidth}{\tablecaption@size{\tablename@size #1}\tablename@skip #2\par}
    \vskip\tablecaption@skip}

%%%%%%%%%%%%%%% TABNOTE
\newcounter{tabnote}
\def\thetabnote{\@fnsymbol\c@tabnote}

\def\tabnotemark@fmt#1{\hbox{$^{#1}$}}
%\def\tabnotemark@fmt#1{\@textsuperscript{#1}}

\def\tabnoteref@sep{,}

% Isvedame zymes
\def\tabnoteref{\@ifnextchar[{\@tempswatrue\@tabnoteref}{\@tempswafalse\@tabnoteref[]}}

\def\@tabnoteref[#1]#2{%
  \if@tempswa% []
    \tabnotemark@fmt{#1}%
  \else%
    \let\@tempa\@empty%
    \tabnotemark@fmt{\@for\@tempb:=#2\do{%
       \@tempa\let\@tempa\tabnoteref@sep%
       \edef\@tempb{\expandafter\@firstofone\@tempb\@empty}%
       \myref{\the\c@table @\@tempb tabnote}}}%
   \fi}

% \tabnotemark neturi prasmes:
\def\tabnotemark{\@ifnextchar[{\@tabnotemark}{\@tabnotemark[]}}

\def\@tabnotemark[#1]#2{\@latex@error{Naudokite komanda \string\tabnoteref!}\@ehc}%

% \@tabnotelabel - surisa zyme su label (kaip ir \thankslabel}

\def\@tabnotelabel[#1]#2{%
   \if@tempswa% []
     \protected@edef\@currentlabel{#1}%
   \else% 
     \refstepcounter{tabnote}%
   \fi%
   \mylabel{\the\c@table @#2tabnote}}%


% \tabnotetext isveda zyme ir teksta
\def\first@tabnote{0}

\def\tabnotetext{\@ifnextchar[{\@tempswatrue\@tabnotetext}{\@tempswafalse\@tabnotetext[]}}

\long\def\@tabnotetext[#1]#2#3{%
  \@tabnotelabel[#1]{#2}%
  \if@tempswa\ifx.#1.\notabnotemark\fi\fi%
  \ifcase\first@tabnote\firsttabnote@skip\def\first@tabnote{1}\fi%
  \tabnote@skip%
  \parbox[t]{\tablewidth}{\tabnote@size\tabnote@indent\tabnotemark@fmt{\expandafter\myref{\the\c@table @#2tabnote}}#3\par}%
  \vtop to\z@{ }\par}

%
\def\notabnotemark{\global\let\sv@tabnotemark@fmt\tabnotemark@fmt\gdef\tabnotemark@fmt##1{\global\let\tabnotemark@fmt\sv@tabnotemark@fmt}}

%%%%%%%%%%%%%%% TABLE (calc. tablewidth)
\newbox\saved@tabular

\newenvironment{Table}
               {\let\@makecaption\@maketablecaption\starttable@hook%
                \let\caption\tabcaption\let\endtabular\endTabular%
                \let\@floatboxreset\table@boxreset\@float{table}}
               {\end@float}

\newenvironment{Table*}
               {\let\@makecaption\@maketablecaption\starttable@hook%
                \let\caption\tabcaption\let\endtabular\endTabular%
                \let\@floatboxreset\table@boxreset\@dblfloat{table}}
               {\end@dblfloat}

\let\Tabular\tabular
\def\endTabular{\crcr\egroup\egroup $\egroup%$
  \setbox\saved@tabular\lastbox
  \global\tablewidth\the\wd\saved@tabular
  \the\cap@temp\global\cap@temp={}\nointerlineskip
  \leavevmode\box\saved@tabular\par}

\newtoks\cap@temp \cap@temp={}
\def\tabcaption{\refstepcounter\@captype \@dblarg{\@tabcaption\@captype}}
\long\def\@tabcaption#1[#2]#3{%
%  \par
%  \addcontentsline{\csname ext@#1\endcsname}{#1}%
%    {\protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
  \global\cap@temp={%
    \begingroup
      \@parboxrestore
      \if@minipage\@setminipage\fi
      \normalsize\@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}\par
    \endgroup}%
  }

%%%%%%%%%%%%%%% HLINE in tabular

   \let\savehline\hline
   \def\thline{\noalign{\vskip3pt}\savehline\noalign{\vskip3pt}}%
   \def\fhline{\noalign{\vskip1pt}\savehline\noalign{\vskip7pt}}%
   \def\bhline{\noalign{\vskip3pt}\noalign{\global\arrayrulewidth=1\p@}\savehline\noalign{\global\arrayrulewidth=.5\p@}\noalign{\vskip3pt}}%
   \def\lhline{\noalign{\vskip3pt}\noalign{\global\arrayrulewidth=.3\p@}\savehline\noalign{\global\arrayrulewidth=.5\p@}\noalign{\vskip3pt}}


%%%%%%%%%%%%%%% SIDEWAYSTABLE, SIDEWAYSFIGURE

\def\local@sidewaystable{%
      \let\@makecaption\@maketablecaption%
      \let\@floatboxreset\table@boxreset%
      \@rotfloat{table}}

\def\local@sidewaysfigure{%
      \let\@makecaption\@makefigurecaption%
      \let\@floatboxreset\figure@boxreset%
      \@rotfloat{figure}}

\AtBeginDocument{\let\sidewaysfigure\local@sidewaysfigure\let\sidewaystable\local@sidewaystable}

%%%% end   of vtexflt.sty

  \renewcommand\figurename{Fig.}
  \setattribute{figurename}{size}{\scshape}
  \setattribute{figurename}{skip}{\hskip1em}
  \setattribute{figurecaption}{size}{\footnotesize\itshape}

  \def\fnum@table{\tablename~\thetable}
  \setattribute{starttable}{hook}{\let\hline\thline}
  \setattribute{tablename}{skip}{\endgraf}
  \setattribute{tablecaption}{size}{\centering\footnotesize\itshape}
  \setattribute{tablename}{size}{\scshape}

  \def\thline{\noalign{\vskip3pt}\savehline\noalign{\vskip4pt}}%
  \def\bhline{\noalign{\vskip3pt}\noalign{\global\arrayrulewidth=1\p@}\savehline\noalign{\global\arrayrulewidth=.5\p@}\noalign{\vskip4pt}}%
  \def\fhline{\noalign{\vskip1pt}\savehline\noalign{\vskip4pt}}%

  \setattribute{firsttabnote}{skip}{\relax}

% ******************** HEADINGS
%
% Normalus RH
\def\ps@headings{%
      \let\@oddfoot\@empty%
      \let\@evenfoot\@empty%
      \def\@evenhead{\runninghead@size\rlap{\pagenumber@size\thepage}\hfill\leftmark\hfill}%
      \def\@oddhead{\runninghead@size\hfill\rightmark\hfill\llap{\pagenumber@size\thepage}}%
}%
% Tuscias RH
\def\ps@empty{\let\@mkboth\@gobbletwo%
     \let\@oddhead\@empty%
     \let\@oddfoot\@empty%
     \let\@evenhead\@oddhead%
     \let\@evenfoot\@oddfoot}
%
% RH  su puslapio numeriu apacioje
\def\ps@plain{\let\@mkboth\@gobbletwo%
     \let\@oddhead\@empty%
     \let\@oddfoot\@empty%
     \let\@evenhead\@oddhead%
     \def\@oddfoot{\hfill\pagenumber@size\thepage\hfill}%
     \let\@evenfoot\@oddfoot}

% First page RH
\def\ps@copyright{\let\@mkboth\@gobbletwo%
  \def\@evenhead{\parbox[t]{\textwidth}{\copyright@size\copyright@text}}%
  \let\@oddhead\@evenhead%
  \def\@oddfoot{\hfill\pagenumber@size\thepage\hfill}%
  \let\@evenfoot\@oddfoot}


% ************************ FOOTNOTE
\def\footnoterule{\kern-5\p@ \hrule width 5.5pc height .2pt depth0pt
  \kern 4.8\p@}                 % The \hrule is hairline

\newcommand\@makefntext[1]{%
    \parindent12pt\@makefnmark #1}
\def\@makefnmark{\@textsuperscript{\normalfont\@thefnmark}}%
%
\def\freefootnotetext[#1]{%
       \begingroup\unrestored@protected@xdef\@thefnmark{#1}\endgroup\@footnotetext}

%
% ******************** Sectioning commands

\newcommand\section{\@startsection {section}{1}{\parindent}%
                                   {\medskipamount}%
                                   {-10pt}%
                                   {\normalsize\upshape\bfseries\mathversion{bold}}}

\newcommand\subsection{\@startsection {subsection}{2}{\parindent}%
                                   {\medskipamount}%
                                   {-10pt}%
                                   {\subsection@shape}}
\def\subsection@shape{\normalsize\itshape}
\def\subsection@prefix{\upshape}

\newcommand\subsubsection{\@startsection{subsubsection}{3}{\parindent}%
                                     {\medskipamount}%
                                     {-10pt}%
                                     {\normalsize\itshape}}
\def\subsubsection@prefix{\upshape}


\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {\smallskipamount}%
                                    {-1em}%
                                    {\normalsize\itshape}}

\newcommand\subparagraph{\@startsection{subparagraph}{5}{\z@}%
                                       {0.1pt}%
                                       {-1em}%
                                       {\normalsize\itshape}}

% Section for the BIBLIOGRAPHY, APPENDIX
\newcommand\specialsection{\@startsection {section}{1}{\z@}%
                                   {\medskipamount}%
                                   {\smallskipamount}%
                                   {\normalsize\centering\MakeUppercase}}

% Format for the counter:
  \def\section@numbersep{.}
  \def\subsection@numbersep{.}
  \def\subsubsection@numbersep{.}
  \def\paragraph@numbersep{.}
  \def\subparagraph@numbersep{.}

  \def\@seccntformat#1{{\csname #1@prefix\endcsname\csname the#1\endcsname\csname#1@numbersep\endcsname\enspace}}

% *********************** APPENDIX

% appendix mess up everything in this layout:


  \newcommand\appendix{\par
    \def\subsection@shape{\normalsize\upshape\bfseries\mathversion{bold}}
    \let\section\specialsection
    \def\section{\@ifnextchar*{\@appsectionstar}{\@appsectionnostar}}%
    \def\section@prefix{\appendixname\ }%
    \def\section@numbersep{:}%
    \setcounter{section}{0}%
    \setcounter{subsection}{0}%
    \gdef\thesection{\@Alph\c@section}}

% \section*{Appendix} -> APPENDIX
  \def\@appsectionstar*#1{%
     \specialsection*{#1}%
     \setcounter{section}{1}%
}

% \section{}      -> APPENDIX A
% \section{proof} -> APPENDIX A: PROOF
  \def\@appsectionnostar#1{%
     \ifx.#1.% 
       \write@toc@ignorenumberline%
       \def\section@numbersep{}%
       \specialsection[\appendixname\ \thesection]{}%
       \write@toc@restorenumberline%
     \else%
       \write@toc@ignorenumberline%
       \def\section@numbersep{:}%
       \specialsection{#1}%
       \write@toc@restorenumberline%
     \fi}


\def\endappendix{}

%****************************** FRONTMATTER *
% FRONT MATTER FORMATTING PARAMETERS

\setattribute{frontmatter} {style} {\centering}
\setattribute{abstract}    {style} {\normaltext}

% FRONT MATTER SKIPS
\setattribute{title}       {skip} {12\p@}
\setattribute{authors}     {skip} {12pt}
\setattribute{affiliation} {skip} {6\p@}
\setattribute{dedicated}   {skip} {12\p@}
\setattribute{address}     {skip} {6\p@ plus 1\p@ minus 1\p@}
\setattribute{abstract}    {skip} {10\p@}
\setattribute{frontmatter} {cmd}  {\vskip20\p@ plus 3\p@ minus 3\p@}

% FRONT MATTER DIMENSIONS
\setattribute{abstract} {width}  {23pc}
\setattribute{abstract} {indent} {16pt} %2em

% FRONT MATTER FONTS 
\setattribute{dochead}    {size} {\Large\bfseries\spaceskip.5em}
\setattribute{dochead}    {case} {\MakeUppercase}
\setattribute{title}      {size} {\normalsize\bfseries\mathversion{bold}\spaceskip.5em}
\setattribute{title}      {case} {\MakeUppercase}
\setattribute{author}     {size} {\normalsize\scshape}
\setattribute{affiliation}{size} {\normalsize\itshape}
\setattribute{dedicated}  {size} {\normalsize\itshape}
\setattribute{abstract}   {size} {\footnotesize\upshape}
\setattribute{address}    {size} {\scriptsize\scshape}
\setattribute{ead}        {size} {\upshape}
\setattribute{history}    {size} {}

\setattribute{runninghead}{size} {\footnotesize}
\setattribute{runninghead}{case} {\MakeUppercase}
\setattribute{pagenumber} {size} {\small}
\setattribute{copyright}  {size} {\fontsize{6}{7}\selectfont\raggedright}

% TEXT, etc.
\setattribute{copyright}  {text} {\url@fmt{}{\itshape}{\@journal}{\journal@url}\\
                                  \the\@pubyear, \volume@label\ \@volume, \issue@label\ \@issue, \@pagerange\\%
                                  \url@fmt{DOI: }{}{\@doi}{\doi@base\@doi}\\
                                  \@ifundefined{@referstodoi}{}{\url@fmt{\referstodoi@prefix}{}{\@referstodoi}{\doi@base\@referstodoi}\\}
                                  \copyright@owner}
\setattribute{copyright} {owner} {$\copyright$~\url@fmt{}{}{Institute of Mathematical Statistics}{\owner@url}, \the\@pubyear}
\setattribute{referstodoi}{prefix}{Main article DOI: }
\setattribute{author}   {prefix} {By~}


\setattribute{arxiv}      {text} {This is an electronic reprint of the original article published by the
                                  \url@fmt{}{}{Institute of Mathematical Statistics}{\owner@url}\ in 
                                  \url@fmt{}{\itshape}{\@journal}{\journal@url},
                                  \url@fmt{}{}{\the\@pubyear, \volume@label\ \@volume, \issue@label\ \@issue, \@pagerange}{\doi@base\@doi}.
                                  This reprint differs from the original in pagination and typographic detail.}
\setattribute{arxiv}     {style} {\small\raggedright}


% HISTORY
\setattribute{history}  {prefix}  {}
\setattribute{history}  {postfix} {.}
\setattribute{received} {prefix}  {Received~}
\setattribute{received} {postfix} {}
\setattribute{revised}  {prefix}  {; revised~}
\setattribute{revised}  {postfix} {}
\setattribute{accepted} {prefix}  {; accepted~}
\setattribute{accepted} {postfix} {}

% COUNTERS, ETC
\newcounter{author}
\newcounter{address}
\newdimen\sv@mathsurround
%
% RUNNING HEAD
\def\runtitle#1{\gdef\@runtitle{\runninghead@case{#1}}}                      \def\@runtitle{}
\def\runauthor#1{{\def\etal{et al.}\gdef\@runauthor{\runninghead@case{#1}}}} \def\@runauthor{}

% DOI
\def\doi#1{\gdef\@doi{#1}}
\def\@doi{}
\def\doi@base{http://dx.doi.org/}

% referstodoi - nurodo i pagrindinio dokumento DOI 
% (is discusion i pagrindini straipsni)
\def\referstodoi#1{\gdef\@referstodoi{#1}}

% \relateddoi - nurodo isvestinius straipsnius
% (is pagrindinio i isvestinius)

\def\relateddois{\@ifnextchar[{\@tempswatrue\@relateddois}{\@tempswafalse\@relateddois[]}}

\def\@relateddois[#1]#2#3{%
  \@thankslabel[#1]{#2}%
  \bgroup%
     \def\relateddoi##1##2{\protect\safe@href{\doi@base##2}{##2}}%
     \protected@xdef\@thanks{\@thanks\protect\thanks@thefnmark{#2thanks}%
     \protect\@footnotetext{\protect\raggedright\thanks@hypertarget{#2}#3}}%
  \egroup%
}%


% \corref
\def\corref#1{}

% \docsubty
\def\docsubty#1{}

%%%%% FRONTMATTER STUFF

\newbox\fm@box
\newdimen\fm@size

\let\hy@frontmatter\relax
\let\hy@endfrontmatter\relax

\def\frontmatter{%
  \global\c@author\z@
  \global\c@address\z@
  \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
%
  \def\pdftitle##1{\write@pdfinfo{\user@hy@title}{##1}}
  \def\pdfauthor##1{\write@pdfinfo{\user@hy@author}{##1}}
  \def\pdfsubject##1{\write@pdfinfo{\user@hy@subject}{##1}}
  \def\pdfkeywords##1{\write@pdfinfo{\user@hy@keywords}{##1}}
%
  \sv@mathsurround\mathsurround \m@th
  \parindent\z@
  \hy@frontmatter
  \open@fm \ignorespaces}
%
%
% ENDFRONTMATTER
%
\def\endfrontmatter{%
  \global\@topnum\z@
  \set@pagerange
  \markboth{\@runauthor}{\@runtitle}%
  \thispagestyle{copyright}%
%
  \history@fmt
  \@thanks
  \keyword@fmt
  \arxivnotice@fmt
  \abstract@fmt
%
  \close@fm
%
  \write@pdfinfo{\hy@author}{\the\authors@list}
  \write@pdfinfo{\hy@subject}{\@journal, \the\@pubyear, \volume@label \@volume, \issue@label \@issue, \@pagerange}
  \write@pdfinfo{\hy@keywords}{\the\keywords@list}
%
  \hy@endfrontmatter
% 
 \global\mathsurround\sv@mathsurround
  \global\c@footnote\z@
  \global\let\@thanks\@empty
  \let\title\relax       \let\author\relax
  \let\address\relax
  \let\frontmatter\relax \let\endfrontmatter\relax
  \let\@maketitle\relax  \let\@@maketitle\relax
%
  \aftergroup\frontmatter@cmd
  }
%
%
\let\maketitle\relax
%
\newdimen\t@xtheight
\def\init@settings{
\splittopskip=\topskip \splitmaxdepth=\maxdepth
\t@xtheight\textheight \advance\t@xtheight-\splittopskip}

\def\no@harm{\let\thanks=\@gobble\let\thanksref=\@gobble\let~\space\def\ead[##1]##2{}\let\\=\@empty \def\protect{\noexpand\protect\noexpand}}

\def\open@fm{%
  \global\setbox\fm@box=\vbox\bgroup
  \hsize=\textwidth
  \frontmatter@style}

\def\close@fm{%
  \par \egroup
  \fm@size=\dp\fm@box \advance\fm@size by \ht\fm@box
  \@whiledim\fm@size>\t@xtheight \do{%
    \global\setbox\@tempboxa=\vsplit\fm@box to \t@xtheight
    \unvbox\@tempboxa %\newpage
    \fm@size=\dp\fm@box \advance\fm@size by \ht\fm@box}
  \if@twocolumn
    \emergencystretch=1pc \twocolumn[\unvbox\fm@box]
  \else
    \unvbox\fm@box
  \fi}

% DOCHEAD
\def\dochead#1{%
  \bgroup
    \dochead@size
    \leavevmode\vphantom{\strut}\dochead@case{#1}\par
  \egroup
  \setattribute{title}{skip}{8\p@}}

% TITLE
\def\title#1{%
  \vglue\title@skip%
  \bgroup%
    \no@harm%
    \let\protect\relax%
    \pdfstringdef\@argi{#1}%
  \egroup%
  \write@pdfinfo{\hy@title}{\@argi}%
  \bgroup%
    \def\\{\break}
    \title@size\title@case{#1}\par%
  \egroup}

%
% AUG - author block
\newenvironment{aug}{}{}

% AUTHOR
%
\def\and{\unskip~and~}

\def\author{\@ifnextchar[{\author@fmt}{\author@fmt[]}}

\def\author@fmt[#1]#2{%
  \stepcounter{author}%
  \ifnum\c@author=1 \author@fmt@init\fi%
  \bgroup% 
     \def\degs##1{##1}\def\fnms##1{##1}\def\inits##1{##1}\def\snm##1{##1}\def\roles##1{##1}%
     \ifnum\c@author=1\addto@authors@list{#2}\else\addto@authors@list{, #2}\fi%
     #2%
  \egroup}

\def\author@fmt@init{%
   \vskip\authors@skip%
   \leavevmode\author@size\author@prefix}

% AFFILIATION
\newcounter{aoasthanksaff}
%\def\@aoas{AOAS}
\def\affiliation#1{%
  \vskip\affiliation@skip
  \bgroup
    %\ifx\@jid\@aoas
      \def\thanksmark{\@ifnextchar[{\@tempswatrue\@thanksmark}{\@tempswatrue\refstepcounter{aoasthanksaff}\@thanksmark[\fnsymbol{aoasthanksaff}]}}%
    %\fi
    \affiliation@size #1\par
  \egroup}


% DEDICATED
\def\dedicated#1{%
  \vskip\dedicated@skip
  \bgroup
    \dedicated@size #1\par
  \egroup}



% ADDRESS
% adresai isvedami straipsnio gale su \printadresses komanda
% adresu tekstus sukisam i token'us
\def\email#1{{e-mail:\ #1}}
%
\def\address{\@ifnextchar[{\address@fmt}{\address@fmt[default]}}
%
% \def\address@ref@A=1
% \newtoks\@address@1
% \@address@1={#2}
%
\def\address@fmt[#1]#2{%
  \stepcounter{address}%
  \expandafter\protected@xdef\csname address@ref@#1\endcsname{\the\c@address}%
  \expandafter\newtoks\csname @address@\the\c@address\endcsname
  \expandafter\expandafter\expandafter\global\csname @address@\the\c@address\endcsname={#2}}
%
% isvedame viena adresa simboliniu vardu:
% printaddress{A}
\def\printaddress#1{%
\bgroup
\xdef\@tmp{\csname address@ref@#1\endcsname}%
\@ifundefined{address@ref@#1}{\@latex@error{Klaida: nera adreso su zyme '#1'!}{??}}{
  \address@size
  \begin{tabular}[t]{@{}l@{}}
  \expandafter\expandafter\expandafter\the\csname @address@\@tmp\endcsname
  \end{tabular}}
\egroup
}

\let\safe@phantomsection\relax

% isvedame viena pagal numeri:
\def\printaddressnum#1{%
\xdef\@tmp{#1}%
\bgroup
\@ifundefined{@address@#1}{\@latex@error{Klaida: nera adreso su numeriu '#1'!}{??}}{
  \address@size
  \ifnum#1=1%
    \write@toc@ignorecontentsline%       
    \safe@phantomsection%
    \addcontentsline{toc}{section}{Author's addresses}%
    \write@toc@restorecontentsline%      
  \fi%
  \begin{tabular}[t]{@{}l@{}}
  \expandafter\expandafter\expandafter\the\csname @address@\@tmp\endcsname
  \end{tabular}}
\egroup
}


% isvedame visus adresus:

\def\address@par{\par\vskip3pt}

\def\printaddresses{%
\vskip\address@skip%
%\addcontentsline{toc}{section}{Author's addresses}%
\def\last@right@glue{\par}%
\ifodd\c@address   \def\last@right@glue{\hfill\hbox{}\address@par} \fi%
\ifnum\c@address=1 \def\last@right@glue{\address@par}\fi%
\@tempcnta=0%
\bgroup\parindent\z@
\@whilenum{\@tempcnta<\c@address}%
  \do{%
    \advance\@tempcnta\@ne
    \ifodd\@tempcnta \def\left@glue{} \def\right@glue{} %     nelyginis
       \else \def\left@glue{\hfill} \def\right@glue{\address@par}\fi % lyginis
    \ifnum\@tempcnta=\c@address \let\left@glue\hfill \let\right@glue\last@right@glue\fi %paskutinis narys
    \left@glue\expandafter\printaddressnum{\the\@tempcnta}\right@glue%
     }
\egroup
}

% ABSTRACT

\newbox\abstract@box

\def\abstract{%
     \global\setbox\abstract@box=\vbox\bgroup%
     \@tempdima\textwidth%
     \advance\@tempdima by-\abstract@width%
     \divide\@tempdima by2%
     \abstract@style%
     \leftskip\@tempdima\rightskip\@tempdima%
     \abstract@size%
     \parindent\sv@parindent%
     \noindent\hskip\abstract@indent\abstractname\ignorespaces}
\def\endabstract{\par\egroup}

\def\abstract@fmt{%
  \ifvoid\abstract@box\else
    \vskip\abstract@skip%
    \unvbox\abstract@box
  \fi}


% HISTORY: received, revised, accepted

\def\history@exist{0}

\def\received#1{\def\@tempa{#1}\ifx\@tempa\@empty\else\gdef\@received{#1}\gdef\history@exist{1}\fi}
  \def\@received{\@nil}
\def\revised#1{\def\@tempa{#1}\ifx\@tempa\@empty\else\gdef\@revised{#1}\gdef\history@exist{1}\fi}
  \def\@revised{\@nil}
\def\accepted#1{\def\@tempa{#1}\ifx\@tempa\@empty\else\gdef\@accepted{#1}\gdef\history@exist{1}\fi}
  \def\@accepted{\@nil}

\def\empty@data{\@nil}
\def\history@fmt{%
  \ifcase\history@exist\else%
  \bgroup
  \protected@xdef\@thefnmark{}%
  \my@footnotetext{%
    \history@size
    \leavevmode
    \history@prefix%
    \ifx\@received\empty@data \else%
      \received@prefix\@received \received@postfix%
    \fi%
    \ifx\@revised\empty@data \else%
      \revised@prefix\@revised \revised@postfix%
    \fi%
    \ifx\@accepted\empty@data \else%
      \accepted@prefix\@accepted \accepted@postfix%
    \fi%
  \history@postfix}
  \egroup
\fi
}

\def\sday#1{#1}
\def\smonth#1{\@ifundefined{month@item@#1}%
     {\@latex@error{Nera tokio menesio, kurio numeris #1!}{??}}%
     {\getitemvalue{month}{#1}}%
  }%
\def\syear#1{#1}
\setvaluelist{month}{January,February,March,April,May,June,July,August,September,October,November,December}


% KEYWORDS

\newtoks\keyword@toks

\setattribute{keyword}{AMS}{\ims@href{\msc@url}{AMS 2000 subject classifications.}}
\setattribute{keyword}{KWD}{Key words and phrases.}
% raktas=class
\setattribute[default]{keyword}{class}{KWD}

%  \keyword@class-> KWD
%  \keyword@KWD  -> AMS 2000... 

% keywords surenkame i token'a
\newpseudoenvironment{keyword}{\gdef\keyword@exist{1}\get@keyword@toks}{}
\def\get@keyword@toks#1\end{\keyword@toks=\expandafter{\the\keyword@toks\keyword#1\endkeyword}\@gobble}

\def\keyword@exist{0}
% \newpseudoenvironment ideda \relax - suvalgome ji
\gdef\keyword#1{\@ifnextchar[{\@keyword}{\@keyword[class=KWD]}}

\gdef\@keyword[#1]{%
  \do@option@list{keyword}{#1}%
  \let\kwd@sep\relax
  \newcounter{kwd\keyword@class}%
  \@ifundefined{kwd\keyword@class num}{\global\@namedef{kwd\keyword@class num}{0}}{}%
  \pre@kwd\hbox{\textit{\csname keyword@\keyword@class\endcsname}\ }%
}

\gdef\endkeyword{%
  \protected@write\@auxout{}{\string\global\string\@namedef{kwd\keyword@class num}{\the\value{kwd\keyword@class}}}%
  \gdef\pre@kwd{\par\leavevmode}}
\let\pre@kwd\relax

% \kwd[; ]{foo}
\def\sep{\unskip\string, }%

\DeclareRobustCommand*\kwd{\stepcounter{kwd\keyword@class}\@ifnextchar[\@kwd{\@kwd[\kwd@sep]}}%
\def\@kwd[#1]#2{\unskip#1{#2}\addto@keywords@list{#2, }\let\kwd@sep\sep%
  \ifnum\value{kwd\keyword@class}=\@nameuse{kwd\keyword@class num}\relax
    \def\check@last@kwd@dot{\@ifnextchar.{}{\unskip.}}%
  \else
    \let\check@last@kwd@dot\relax%
  \fi\check@last@kwd@dot%
}%  


\def\keyword@fmt{\ifcase\keyword@exist\else\bgroup\protected@xdef\@thefnmark{}\my@footnotetext{\the\keyword@toks}\egroup\fi}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% supplement
%
%\begin{supplement}[id=suppA]
%  \sname{Supplement A}
%  \stitle{}
%  \slink[doi]{}
%  \slink[url]{http://lib.stat.cmu.edu/aoas/???/???}
%  \sdescription{}
%\end{supplement}
% \thesuppdoi{suppA}
% \ref{suppA}
% \hyperlink{suppA}{text}

\def\supplement@name{Supplementary Material}

\def\sname#1{\def\@sname{#1}\def\@currentlabel{#1}}
\def\stitle#1{\def\@stitle{#1}}
\def\sdatatype#1{\def\@sdatatype{#1}}
\def\sfilename#1{\def\@sfilename{}}
\define@key{slink}{text}{\def\sdoi@text{#1}}
\define@key{slink}{doi}[true]{\def\slink@type{doi}}
\define@key{slink}{url}[true]{\def\slink@type{url}}
\def\slink[#1]#2{\setkeys{slink}{#1}\@namedef{slink@\slink@type}{#2}}
\def\sdescription#1{\def\@sdescription{#1}}

\def\suppsection@fmt{\specialsection*{\supplement@name}}

\def\slink@doi@fmt{%
  \url@fmt{DOI: }{}{\slink@doi}{\doi@base\slink@doi}%
  \@ifundefined{supp@label}{}{%
    \expandafter\xdef\csname\supp@label @doi\endcsname{\slink@doi}%
    \expandafter\xdef\csname\supp@label @url\endcsname{\doi@base\slink@doi}}%
}

\def\thesuppdoi#1{\@ifundefined{#1@doi}%
  {\@latex@error{Neapibreztas supplement id=#1}{??}}%
  {\def\@tempx{\csname #1@doi\endcsname}%
   \@ifundefined{#1@url}{\def\@tempy{\doi@base\csname #1@url\endcsname}}{\def\@tempy{\csname #1@url\endcsname}}%
   \url@fmt{DOI: }{}{\@tempx}{\@tempy}}%
}

\define@key{supplement}{id}{\def\supp@label{#1}}

\long\def\supplement{\@ifnextchar[{\@supplement}{\@supplement[]}}

\long\def\@supplement[#1]{%
  \suppsection@fmt%
  \global\let\suppsection@fmt\smallskip%
  \setkeys{supplement}{#1}%
}

\def\endsupplement{%
  \@ifundefined{@sname}{}{\@ifundefined{supp@label}{}{\hypertarget{\supp@label}{}\label{\supp@label}}\textbf{\@sname: }}%
  \@ifundefined{@stitle}{}{\textbf{\@stitle\ }}%
%
  (%
  \@ifundefined{slink@doi}{}{\slink@doi@fmt}%
  \@ifundefined{@sdatatype}{}{; \@sdatatype}%
  \@ifundefined{@sfilename}{}{\@sfilename}%
  ).
  \@ifundefined{@sdescription}{}{\@sdescription}%
  \par}


%%%%%%%%%  special notice of arXiv:

\def\arxivnotice@fmt{%
  \bgroup
    \protected@xdef\@thefnmark{}
    \my@footnotetext{\vskip-3\p@\noindent\fbox{\parbox{\textwidth}{\arxiv@style\arxiv@text}}}
  \egroup}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \tableofcontents

%\RequirePackage{imstoc}

%%%%%%%%%% contents of imspdf  [2008/12/11]


%%%%%%%%%% hooks:
\setattribute{toc}    {size}{\footnotesize}
\setattribute{tocname}{size}{\footnotesize}

\setattribute{toc}     {skip}{\relax}
\setattribute{tocname} {skip}{\relax}
\setattribute{aftertoc}{skip}{\relax}
\setattribute{toc}{heading}{\nocontentsline\section*{\contentsname}} 

\setattribute{tocline}{dot} {.} 
\setattribute{tocline}{skip}{\enskip} 

%%%%%%%%%% papildomos komandos pagrindinem klasem:

% ignoruoja numberline (appendix)
\def\ignorenumberline{\let\sv@numberline\numberline\let\numberline\@gobble}
\def\restorenumberline{\let\numberline\sv@numberline}

\def\write@toc@ignorenumberline {\addtocontents{toc}{\protect\ignorenumberline}}
\def\write@toc@restorenumberline{\addtocontents{toc}{\protect\restorenumberline}}

% neitraukia i turini (printaddress)

\def\ignorecontentsline{\let\sv@contentsline\contentsline\let\contentsline\@gobblefour}
\def\restorecontentsline{\let\contentsline\sv@contentsline}

% norim kad bookmark'uose butu, TOC - ne
\def\write@toc@ignorecontentsline {\addtocontents{toc}{\protect\ignorecontentsline}}
\def\write@toc@restorecontentsline{\addtocontents{toc}{\protect\restorecontentsline}} 

%%%%%%%%%% parametrai


\newcommand\@pnumwidth{2.00em}%1.55em
\newcommand\@tocrmarg{2.55em plus 1fil}
\newcommand\@dotsep{4.2}
%\newcommand\@dotsep{10000}

\def\contentsname{Contents}


\define@key{tableofcontents}{level}{%
  \@tempcnta=#1\relax%
  \ifcase\@tempcnta
     \@latex@error{TOC lygio =0 nera!}{}
  \or% 1
  \or% 2
     \let\l@subsection\lx@subsection
  \or% 3
     \let\l@subsection\lx@subsection
     \let\l@subsubsection\lx@subsubsection
  \or% 4
    \@latex@error{TOC lygio =4 nera!}{}
  \fi}

% sekciju numerius lyginame kaireje
\define@key{tableofcontents}{alignleft}[true]{\let\numberline\left@numberline}

% sekciju numerius lyginame desneje
\define@key{tableofcontents}{alignright}[true]{\let\numberline\right@numberline}

\def\tableofcontents{\@ifnextchar[{\@tableofcontents}{\@tableofcontents[]}}

\def\@tableofcontents[#1]{%
    \setkeys{tableofcontents}{#1}%
    \tocname@skip%
    \bgroup
      \tocname@size%
      \toc@heading%
      \toc@size%
      \toc@skip%
      \@starttoc{toc}%
    \egroup
    \aftertoc@skip
    }

\newdimen\cw@section
\newdimen\cw@subsection
\newdimen\cw@subsubsection

\newdimen\wd@subsubsection

\def\maxcsec#1{%
  \setbox\@tempboxa\hbox{\toc@size#1\tocline@dot\tocline@skip}%
  \global\cw@section\wd\@tempboxa}

\def\maxcsubsec#1{%
  \setbox\@tempboxa\hbox{\toc@size#1\tocline@dot\tocline@skip}%
  \global\cw@subsection\wd\@tempboxa%
}


\def\maxcsubsubsec#1{%
  \setbox\@tempboxa\hbox{\toc@size#1\tocline@dot\tocline@skip}%
  \global\cw@subsubsection\wd\@tempboxa%
  \global\wd@subsubsection\cw@section%
  \advance\wd@subsubsection by\cw@subsection%
}


\AtBeginDocument{\maxcsec{1}\maxcsubsec{1.1}\maxcsubsubsec{1.1.1}}
 
%% keiciame [3] i [2]
\newcommand*\l@section[2]{\@dottedtocline{1}{0em}{\cw@section}{#1}{#2}}
\newcommand*\l@heading[2]{\@dottedtocline{1}{0em}{\cw@section}{#1}{#2}}

\newcommand*\l@subsection[2]{}
\newcommand*\l@subsubsection[2]{}
\newcommand*\l@paragraph[2]{}
\newcommand*\l@subparagraph[2]{}

\newcommand*\lx@subsection[2]{\@dottedtocline{1}{\cw@section}{\cw@subsection}{#1}{#2}}
\newcommand*\lx@subsubsection[2]{\@dottedtocline{1}{\wd@subsubsection}{\cw@subsubsection}{#1}{#2}}

\def\right@numberline#1{\ifx.#1.\else\hb@xt@\@tempdima{\hfil{#1}\tocline@dot\tocline@skip}\fi}

\def\left@numberline#1{\ifx.#1.\else\hb@xt@\@tempdima{{#1}\tocline@dot\tocline@skip}\fi}

\let\numberline\right@numberline

%\def\tocnumber@fmt#1#2{\hyperlink{#2}{#1}}

%\def\notocnumbers{\let\tocnumber@fmt\@gobbletwo}


%%%%%%%%%% end contents of imspdf 


%%

\setattribute{toc}{heading}{\nocontentsline\specialsection*{\contentsname}} 
\setattribute{aftertoc}{skip}{\medbreak}

% appendix komandoje:
\def\ignorenumberline{\let\sv@numberline\numberline\let\numberline\@gobble}
\def\restorenumberline{\let\numberline\sv@numberline}

% printaddress

\def\ignorecontentsline{\let\sv@contentsline\contentsline\let\contentsline\@gobblefour}
\def\restorecontentsline{\let\contentsline\sv@contentsline}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interaction with hyperref

\def\pdfstringdef#1#2{}

\def\test@hyperref{\@ifundefined{Hy@SetCatcodes}{}{\imsart@hyperref@settings}}

\AtBeginDocument{\test@hyperref}

\def\imsart@hyperref@settings{%
% activate href
   \let\ims@href\href%
   \let\safe@href\href
   \let\safe@phantomsection\phantomsection
%
   \def\thankref@hyperlink##1{\edef\@tempx{##1thanks}\hyperlink{##1}{\saferef{\@tempx}}}
   \def\thanks@hypertarget##1{\smash{\raise\baselineskip\hbox{\protect\hypertarget{##1}{}}}}
   \def\thanksmark@link##1##2{\hyperlink{##1}{##2}}
%
% put document info
   \def\write@pdfinfo##1##2{\protected@write\@auxout{\no@harm}{\string\gdef\string##1{##2}}}
   \@ifundefined{hy@title}{}{\pdfstringdef\@pdftitle{\hy@title}}
   \@ifundefined{hy@author}{}{\pdfstringdef\@pdfauthor{\hy@author}}
   \@ifundefined{hy@subject}{}{\pdfstringdef\@pdfsubject{\hy@subject}}
   \@ifundefined{hy@keywords}{}{\pdfstringdef\@pdfkeywords{\hy@keywords}}
%
   \@ifundefined{user@hy@title}{}{\global\let\@pdftitle\user@hy@title}
   \@ifundefined{user@hy@author}{}{\global\let\@pdfauthor\user@hy@author}
   \@ifundefined{user@hy@subject}{}{\global\let\@pdfsubject\user@hy@subject}
   \@ifundefined{user@hy@keywords}{}{\global\let\@pdfkeywords\user@hy@keywords}
%
   \def\saferef##1{\expandafter\safe@setref\csname r@##1\endcsname\@firstoffive{##1}}
%
% Define colors: 
  \hypersetup{
     linktocpage=true,
     colorlinks=true,
     linkcolor=blue,
     citecolor=blue,
     urlcolor =blue
             }
% MathSciNet:
   \def\MR##1{\href{http://www.ams.org/mathscinet-getitem?mr=##1}{MR##1}}
}

\def\write@pdfinfo#1#2{}

\newtoks\authors@list
\def\addto@authors@list#1{%
  \begingroup%
    \no@harm%
    \xdef\@act{\global\noexpand\authors@list{\the\authors@list#1}}\@act%
  \endgroup}

\newtoks\keywords@list
\def\addto@keywords@list#1{%
  \begingroup%
    \no@harm%
    \xdef\@act{\global\noexpand\keywords@list{\the\keywords@list#1}}\@act%
  \endgroup}

% MathSciNet:
\def\MR#1{MR#1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% msnbibl.sty
\if@MSNbibl

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bibitem

  \def\endbibitem{\set@msnbibl\let\endbibitem\relax}

  \def\set@msnbibl{\protected@write\@auxout{}{\string\global\string\let\string\bibitem\string\msn@bibitem}{}}

  \let\real@bibitem\bibitem

  \long\def\msn@bibitem#1\endbibitem{%
     \real@bibitem#1\endbibitem
    }%


  \def\bbl@MR#1{#1}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get@MR

  \def\get@MR#1R#2 #3\end{\MR{#2}}

  \let\MR@fmt\@firstofone

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Common macros:

  % Setting a "style" for a command:
  % \set@bibl@cmd{bvolume}  == \def\bvolume#1{{\bvolume@style #1}}

  \def\set@bibl@cmd#1{\expandafter\def\csname #1\endcsname##1{{\csname #1@style\endcsname##1}}}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  \def\get@battribute#1{\ifx.#1.\else\def\b@attribute{#1}\fi}

  \let\show@battribute\relax

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bauthor, beditor

  \def\bbl@bauthor#1{{%
    \let\binits\@firstofone%
    \let\bsnm\@firstofone%
    \let\bfnm\@gobble%
    \let\bparticle\@firstofone%
    \let\bsuffix\@firstofone%
    \bauthor@style%
  #1}}

  \def\bbl@beditor#1{{%
    \let\binits\@firstofone%
    \let\bsnm\@firstofone%
    \let\bfnm\@gobble%
    \let\bparticle\@firstofone%
    \let\bsuffix\@firstofone%
    \beditor@style%
  #1}}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bid
  % \bid{MR={},doi={},...}

  \define@key{bid}{mr}{\get@MR MR#1 \end\ignorespaces}
  \define@key{bid}{doi}{}
  \define@key{bid}{pubmed}{}
  \define@key{bid}{pii}{}
  \define@key{bid}{issn}{}
  \define@key{bid}{pmid}{}
  \define@key{bid}{pmcid}{}
  \define@key{bid}{mid}{}

  \def\bbl@bid#1{%
         \setkeys{bid}{#1}%
  }

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% common@pub@types

  \def\common@pub@types{%
    \def\AND{and }%
    \let\betal\@firstofone%
    \set@bibl@cmd{btitle}
    \let\byear\@firstofone%  
    \let\bpages\@firstofone%
    \let\bmisc\@firstofone%
    \let\bnote\@firstofone%
    \let\banumber\@firstofone%
    \def\bmrnumber##1{\get@MR##1 \end\ignorespaces}%
    \let\bid\bbl@bid%
    \let\bptnote\@gobble% 
    \set@bibl@cmd{bvolume}%
    \csname common@pub@types@hook\endcsname%
  }


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% default stiliai

    \setattribute{bauthor}{style}{\scshape}
    \setattribute{beditor}{style}{\scshape}

    \setattribute{bjournal}  {style}{\itshape}
    \setattribute{bbooktitle}{style}{\itshape}
    \setattribute{bseries}   {style}{\itshape}

    \setattribute{bvolume}   {style}{\bfseries}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  barticle

  \def\barticle{\@ifnextchar[{\@barticle}{\@barticle[]}}

  \def\@barticle[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \set@bibl@cmd{bjournal}%
  }

  \def\endbarticle{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bbook

  \def\bbook{\@ifnextchar[{\@bbook}{\@bbook[]}}

  \def\@bbook[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\beditor\bbl@beditor%
    \let\bedition\@firstofone%
    \set@bibl@cmd{bseries}%
    \let\bpublisher\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
    \let\bisbn\@gobble%
  %
    \let\btitle@style\itshape
  }

  \def\endbbook{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bincollection

  \def\bincollection{\@ifnextchar[{\@bincollection}{\@bincollection[]}}

  \def\@bincollection[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\beditor\bbl@beditor%
    \set@bibl@cmd{bbooktitle}%
    \let\bchapter\@firstofone%
    \let\bedition\@firstofone%
    \set@bibl@cmd{bseries}%
    \let\bpublisher\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
    \let\bisbn\@gobble%
  }

  \def\endbincollection{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% binproceedings

  \def\binproceedings{\@ifnextchar[{\@binproceedings}{\@binproceedings[]}}

  \def\@binproceedings[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\beditor\bbl@beditor%
    \set@bibl@cmd{bbooktitle}%
    \set@bibl@cmd{bseries}%
    \let\borganization\@firstofone%
    \let\bpublisher\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
    \let\bisbn\@gobble%
  }

  \def\endbinproceedings{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% binbook

  \def\binbook{\@ifnextchar[{\@binbook}{\@binbook[]}}

  \def\@binbook[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\beditor\bbl@beditor%
    \set@bibl@cmd{bbooktitle}%
    \let\bchapter\@firstofone%
    \let\bedition\@firstofone%
    \set@bibl@cmd{bseries}%
    \let\bpublisher\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
    \let\bisbn\@gobble%
  %
  %  \let\btitle@style\itshape
  }

  \def\endbinbook{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bproceedings

  \def\bproceedings{\@ifnextchar[{\@bproceedings}{\@bproceedings[]}}

  \def\@bproceedings[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\beditor\bbl@beditor%
    \let\bedition\@firstofone%
    \set@bibl@cmd{bseries}%
    \let\bpublisher\@firstofone%
    \let\borganization\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
    \let\bisbn\@gobble%
  %
    \let\btitle@style\itshape
  }

  \def\bproceedings{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% btechreport

  \def\btechreport{\@ifnextchar[{\@btechreport}{\@btechreport[]}}

  \def\@btechreport[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\btype\@firstofone%
    \let\bnumber\@firstofone%
    \let\binstitution\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
  }

  \def\endbtechreport{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bmanual

  \def\bmanual{\@ifnextchar[{\@bmanual}{\@bmanual[]}}

  \def\@bmanual[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\bedition\@firstofone%
    \let\borganization\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
    \let\bisbn\@gobble%
  }

  \def\endbmanual{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mastersthesis

  \def\bmastersthesis{\@ifnextchar[{\@bmastersthesis}{\@bmastersthesis[]}}

  \def\@bmastersthesis[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\bschool\@firstofone%
    \let\btype\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
  }

  \def\endbmastersthesis{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% phdthesis

  \def\bphdthesis{\@ifnextchar[{\@bphdthesis}{\@bphdthesis[]}}

  \def\@bphdthesis[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\bschool\@firstofone%
    \let\btype\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
  }

  \def\endbphdthesis{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bbooklet

  \def\bbooklet{\@ifnextchar[{\@bbooklet}{\@bbooklet[]}}

  \def\@bbooklet[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\bhowpublished\@firstofone%
    \let\baddress\@firstofone%
    \let\blocation\@firstofone%
  }

  \def\endbbooklet{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bunpublished

  \def\bunpublished{\@ifnextchar[{\@bunpublished}{\@bunpublished[]}}

  \def\@bunpublished[#1]{%
    \get@battribute{#1}%
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
  }

  \def\endbunpublished{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bmisc

  \def\bmisc{\@ifnextchar[{\@bmisc}{\@bmisc[]}}

  \def\@bmisc[#1]{%
    \get@battribute{unstr}%  
    \common@pub@types%
    \let\bauthor\bbl@bauthor%
    \let\bhowpublished\@firstofone%
    \let\borganization\@firstofone%
  }

  \def\endbmisc{\show@battribute\aftergroup\ignorespaces}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \bptok

  \def\bptok#1{\ignorespaces}

\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \OrigBibText ... \endOrigBibText
\ifx\OrigBibText\undefined \long\def\OrigBibText#1\endOrigBibText{}\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vtexurl.sty

%%%%%%%%%%%%%%%%%%% \doiurl

%% \doiurl{http://dx.doi.org/00000.0000}
%% \doiurl{00000.0000}
%% \doiurl{http://\\dx.doi.org/00000.0000}
%% \doiurl{00000.\\0000}

\@ifundefined{doi@base}{\def\doi@base{http://dx.doi.org/}}{}

\DeclareRobustCommand\doiurl[1]{{%
 \def\\{}%
 \check@doiurl@prefix#1http://\end%
 \check@doiurl@break#1\\\end\relax%
 \if@doiurlbreak%
    \@ifundefined{doi@url}{%
       \href{\doi@base\@tempx\@tempy}{\csname doi@size\endcsname\@tempx}%
       \break%
       \href{\doi@base\@tempx\@tempy}{\csname doi@size\endcsname\@tempy}}%
                          {%
       \expandafter\href{\@tempx\@tempy}{\csname doi@size\endcsname\@tempx}%
       \break%
       \href{\@tempx\@tempy}{\csname doi@size\endcsname\@tempy}}%
 \else%
    \@ifundefined{doi@url}%
        {\href{\doi@base#1}{\csname doi@size\endcsname #1}}%
        {\href{#1}{\csname doi@size\endcsname #1}}%
 \fi}}


\def\check@doiurl@prefix#1http://#2\end{\ifx.#2.\else\def\doi@url{fullurl}\fi}

\newif\if@doiurlbreak \@doiurlbreakfalse

\def\check@doiurl@break#1\\#2\end{%
    \ifx.#2.\@doiurlbreakfalse%
    \else%
       \@doiurlbreaktrue%
       \def\@tempx{#1}%
       \def\@tempy{#2}%
    \fi}

%%%%%%%%%%%%%%%%%%% arxivurl

%% \arxivurl{http://arxiv.org/abs/math.PR/0603300}
%% \arxivurl{math.PR/0603300}
%% \arxivurl{http://\\arxiv.org/abs/math.PR/0603300}
%% \arxivurl{math.PR/\\0603300}

\@ifundefined{arxiv@base}{\def\arxiv@base{http://arxiv.org/abs/}}{}

\DeclareRobustCommand\arxivurl[1]{{%
 \def\\{}%
 \check@arxivurl@prefix#1http://\end%
 \check@arxivurl@break#1\\\end\relax%
 \if@arxivurlbreak%
    \@ifundefined{arxiv@url}{%
       \href{\arxiv@base\@tempx\@tempy}{\csname arxivurl@size\endcsname\@tempx}%
       \break%
       \href{\arxiv@base\@tempx\@tempy}{\csname arxivurl@size\endcsname\@tempy}}%
                            {%
       \expandafter\href{\@tempx\@tempy}{\csname arxivurl@size\endcsname\@tempx}%
       \break%
       \href{\@tempx\@tempy}{\csname arxivurl@size\endcsname\@tempy}}%
 \else%
    \@ifundefined{arxiv@url}%
        {\href{\arxiv@base#1}{\csname arxivurl@size\endcsname #1}}%
        {\href{#1}{\csname arxivurl@size\endcsname #1}}%
 \fi}}


\def\check@arxivurl@prefix#1http://#2\end{\ifx.#2.\else\def\arxiv@url{fullurl}\fi}

\newif\if@arxivurlbreak \@arxivurlbreakfalse

\def\check@arxivurl@break#1\\#2\end{%
    \ifx.#2.\@arxivurlbreakfalse%
    \else%
       \@arxivurlbreaktrue%
       \def\@tempx{#1}%
       \def\@tempy{#2}%
    \fi}


%%%%%%%%%%%%%%%%%%% surl

%% \surl{http://www.foo.com/something}
%% \surl{ftp://ftp.foo.com/something}
%% \surl{www.foo.com/something}
%% \surl{http://www.foo.com/\\something}
%% \surl{ftp://\\ftp.foo.com/\\something}
%% \surl{www.foo.com/\\something}

\def\spliturl@prefix{http://}

\DeclareRobustCommand\surl[1]{{%
 \def\\{}%
 \check@spliturl@prefix#1://\end%
 \check@spliturl@break#1\\\end\relax%
 \if@spliturlbreak%
    \@ifundefined{spliturl@url}{%
       \href{\spliturl@prefix\@tempx\@tempy}{\csname surl@size\endcsname\@tempx}%
       \break%
       \href{\spliturl@prefix\@tempx\@tempy}{\csname surl@size\endcsname\@tempy}}%
                            {%
       \expandafter\href{\@tempx\@tempy}{\csname surl@size\endcsname\@tempx}%
       \break%
       \href{\@tempx\@tempy}{\csname surl@size\endcsname\@tempy}}%
 \else%
    \@ifundefined{spliturl@url}%
        {\href{\spliturl@prefix#1}{\csname surl@size\endcsname #1}}%
        {\href{#1}{\csname surl@size\endcsname #1}}%
 \fi}}


\def\check@spliturl@prefix#1://#2\end{\ifx.#2.\else\def\spliturl@url{fullurl}\fi}

\newif\if@spliturlbreak \@spliturlbreakfalse

\def\check@spliturl@break#1\\#2\end{%
    \ifx.#2.\@spliturlbreakfalse%
    \else%
       \@spliturlbreaktrue%
       \def\@tempx{#1}%
       \def\@tempy{#2}%
    \fi}


%%%%%%%%%%%%%%%%%%% shref
%% \shref{http://www.foo.com/something}{http://\\www.\\foo.\\com/something}

\DeclareRobustCommand\shref[2]{%
 \def\split@href@prefix{}%
 \fmt@split@href[#1]#2\\\\%
 }


\def\fmt@split@href[#1]#2\\{%
   \ifx.#2.\def\next[##1]#2{}%
   \else \split@href@prefix\href{#1}{\csname surl@size\endcsname#2}%
         \gdef\split@href@prefix{\break}%
         \let\next\fmt@split@href
   \fi\next[#1]}

%%%%%%%%%%%%%%%%%% semail (split email)
%% \semail{mailto@foo.\\something.\\lt}

%%
\def\semail@prefix{mailto:}
\DeclareRobustCommand\semail[1]%
{%
 \bgroup
  \let\\\@empty
  \def\_{\string_}%
  \xdef\@tempa{#1}%
 \egroup%
 \shref{\semail@prefix\@tempa}{#1}%
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************** nereikalingi apibrezimai:
\let\HPROOF\relax
\let\PROOF\relax
\let\CRC\relax
\let\psdraft\relax
\let\psfull\relax
\let\vtexed\@gobble

%**************** INICIALIZATION

\@twosidetrue
\pagenumbering{arabic}
\frenchspacing
\init@settings
\pagestyle{headings}

\endinput
%%
%% End of file `ximsims.cls'.
