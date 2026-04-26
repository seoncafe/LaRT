%% kiseon_memo.cls
%% Document class equivalent of kiseon_memo.tex.
%% Usage: \documentclass[<options>]{kiseon_memo}
%%
%% Options passed through to the underlying 'article' class (e.g. a4paper,
%% 10pt, 11pt, 12pt, english, korean, ...).  The class loads its own fixed
%% layout, so font-size options are accepted but may be overridden by scrextend.
%%
%% Packages automatically loaded:
%%   graphicx, amssymb, amsmath, mathrsfs,
%%   fancyhdr, lastpage, kotex,
%%   hyperref (with hyperfootnotes=false),
%%   scrextend (fontsize=10pt)
%%
%% Font: Times Roman for text, CM for Greek/math symbols (from apjfonts.sty
%%   by Alexey Vikhlinin & Maxim Markevitch, 1996; LaTeX2e portion inlined).
%%   Requires ptmcm virtual fonts (part of psnfss/mathptmx).
%%
%% \title, \author, \date behave exactly as in the original kiseon_memo.tex:
%%   they store values used in the header/footer and the boxed \maketitle.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{kiseon_memo}[2026/04/26 v1.1 Kiseon memo class]

%% ---- Pass all unknown options to 'article' --------------------------------
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass{article}

%% ---- Required packages ----------------------------------------------------
\RequirePackage[T1]{fontenc}
%% Fix \AA broken when Times redefines \r (apjfonts.sty workaround)
\let\rold\r
\def\AA{\rold A}
\let\r\relax

\RequirePackage{mathrsfs}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{graphicx}
\RequirePackage{fancyhdr}
\RequirePackage{kotex}
\RequirePackage[hyperfootnotes=false]{hyperref}
\RequirePackage{lastpage}
\RequirePackage[fontsize=10pt]{scrextend}
\RequirePackage[table]{xcolor}

%% ---- APJ Times + CM math fonts (inlined from apjfonts.sty) ----------------
%% Text fonts: Times Roman, Helvetica, Courier
\renewcommand{\sfdefault}{phv}
\renewcommand{\rmdefault}{ptm}
\renewcommand{\ttdefault}{pcr}

%% Symbol font declarations: Times (ptmcm) for math letters, CM for greek/symbols
\DeclareSymbolFont{cmoperators}   {OT1}{cmr}{m}{n}
\DeclareSymbolFont{cmitoperators} {OT1}{cmr}{m}{it}
\DeclareSymbolFont{operators}     {OT1}{ptm}{m}{n}
\DeclareSymbolFont{letters}       {OML}{ptmcm}{m}{it}
\DeclareSymbolFont{ptmletters}    {OML}{ptmcm}{m}{it}
\DeclareSymbolFont{cmletters}     {OML}{cmm}{m}{it}
\DeclareSymbolFont{symbols}       {OMS}{cmsy}{m}{n}
\DeclareSymbolFont{ptmsymbols}    {U}{psy}{m}{n}
\DeclareSymbolFont{largesymbols}  {OMX}{cmex}{m}{n}
\DeclareSymbolFont{bold}          {OT1}{ptm}{bx}{n}
\DeclareSymbolFont{italic}        {OT1}{ptm}{m}{it}

\SetSymbolFont{cmoperators}   {bold}{OT1}{cmr}{b}{n}
\SetSymbolFont{cmitoperators} {bold}{OT1}{cmr}{bx}{it}
\SetSymbolFont{operators}     {bold}{OT1}{ptm}{bx}{n}
\SetSymbolFont{letters}       {bold}{OML}{ptmCM}{b}{it}
\SetSymbolFont{ptmletters}    {bold}{OML}{ptmCM}{b}{it}
\SetSymbolFont{cmletters}     {bold}{OML}{cmm}{b}{it}
\SetSymbolFont{symbols}       {bold}{OMS}{cmsy}{b}{n}
\SetSymbolFont{largesymbols}  {bold}{OMX}{cmex}{m}{n}
\SetSymbolFont{bold}          {bold}{OT1}{ptm}{bx}{n}
\SetSymbolFont{italic}        {bold}{OT1}{ptm}{b}{it}

\DeclareMathAlphabet{\mathbf}{OT1}{ptm}{bx}{n}
\DeclareMathAlphabet{\mathit}{OT1}{ptm}{m}{it}
\DeclareMathAlphabet{\mathbit}{OT1}{ptm}{bx}{it}

\DeclareMathSymbol{\omicron}{\mathord}{operators}{`\o}

\thinmuskip=2mu
\medmuskip=2.5mu plus 1mu minus 1mu
\thickmuskip=4mu plus 1.5mu minus 1mu

\DeclareMathSymbol{-}{\mathbin}{ptmsymbols}{45}
\DeclareMathSymbol{+}{\mathbin}{ptmsymbols}{43}

%% Lowercase Greek — CM letters
\DeclareMathSymbol{\alpha}      {\mathord}{cmletters}{'013}
\DeclareMathSymbol{\beta}       {\mathord}{cmletters}{'014}
\DeclareMathSymbol{\gamma}      {\mathord}{cmletters}{'015}
\DeclareMathSymbol{\delta}      {\mathord}{cmletters}{'016}
\DeclareMathSymbol{\epsilon}    {\mathord}{cmletters}{'017}
\DeclareMathSymbol{\zeta}       {\mathord}{cmletters}{'020}
\DeclareMathSymbol{\eta}        {\mathord}{cmletters}{'021}
\DeclareMathSymbol{\theta}      {\mathord}{cmletters}{'022}
\DeclareMathSymbol{\iota}       {\mathord}{cmletters}{'023}
\DeclareMathSymbol{\kappa}      {\mathord}{cmletters}{'024}
\DeclareMathSymbol{\lambda}     {\mathord}{cmletters}{'025}
\DeclareMathSymbol{\mu}         {\mathord}{cmletters}{'026}
\DeclareMathSymbol{\nu}         {\mathord}{cmletters}{'027}
\DeclareMathSymbol{\xi}         {\mathord}{cmletters}{'030}
\DeclareMathSymbol{\pi}         {\mathord}{cmletters}{'031}
\DeclareMathSymbol{\rho}        {\mathord}{cmletters}{'032}
\DeclareMathSymbol{\sigma}      {\mathord}{cmletters}{'033}
\DeclareMathSymbol{\tau}        {\mathord}{cmletters}{'034}
\DeclareMathSymbol{\upsilon}    {\mathord}{cmletters}{'035}
\DeclareMathSymbol{\phi}        {\mathord}{cmletters}{'036}
\DeclareMathSymbol{\chi}        {\mathord}{cmletters}{'037}
\DeclareMathSymbol{\psi}        {\mathord}{cmletters}{'040}
\DeclareMathSymbol{\omega}      {\mathord}{cmletters}{'041}
\DeclareMathSymbol{\varepsilon} {\mathord}{cmletters}{'042}
\DeclareMathSymbol{\vartheta}   {\mathord}{cmletters}{'043}
\DeclareMathSymbol{\varomega}   {\mathord}{cmletters}{'044}
\DeclareMathSymbol{\varrho}     {\mathord}{cmletters}{'045}
\DeclareMathSymbol{\varsigma}   {\mathord}{cmletters}{'046}
\DeclareMathSymbol{\varphi}     {\mathord}{cmletters}{'047}
\DeclareMathSymbol{\partial}    {\mathord}{cmletters}{'100}
\DeclareMathSymbol{\ell}        {\mathord}{cmletters}{'140}
\def\hbar{{\raise-0.15truept\hbox{$\mkern-0.5mu\barbar$}\mkern-8.5muh}}

%% Uppercase Greek — CM operators
\DeclareMathSymbol\Gamma   {\mathord}{cmoperators}{"00}
\DeclareMathSymbol\Delta   {\mathord}{cmoperators}{"01}
\DeclareMathSymbol\Theta   {\mathord}{cmoperators}{"02}
\DeclareMathSymbol\Lambda  {\mathord}{cmoperators}{"03}
\DeclareMathSymbol\Xi      {\mathord}{cmoperators}{"04}
\DeclareMathSymbol\Pi      {\mathord}{cmoperators}{"05}
\DeclareMathSymbol\Sigma   {\mathord}{cmoperators}{"06}
\DeclareMathSymbol\Upsilon {\mathord}{cmoperators}{"07}
\DeclareMathSymbol\Phi     {\mathord}{cmoperators}{"08}
\DeclareMathSymbol\Psi     {\mathord}{cmoperators}{"09}
\DeclareMathSymbol\Omega   {\mathord}{cmoperators}{"0A}

\DeclareMathSymbol\varGamma   {\mathord}{cmletters}{"00}
\DeclareMathSymbol\varDelta   {\mathord}{cmletters}{"01}
\DeclareMathSymbol\varTheta   {\mathord}{cmletters}{"02}
\DeclareMathSymbol\varLambda  {\mathord}{cmletters}{"03}
\DeclareMathSymbol\varXi      {\mathord}{cmletters}{"04}
\DeclareMathSymbol\varPi      {\mathord}{cmletters}{"05}
\DeclareMathSymbol\varSigma   {\mathord}{cmletters}{"06}
\DeclareMathSymbol\varUpsilon {\mathord}{cmletters}{"07}
\DeclareMathSymbol\varPhi     {\mathord}{cmletters}{"08}
\DeclareMathSymbol\varPsi     {\mathord}{cmletters}{"09}
\DeclareMathSymbol\varOmega   {\mathord}{cmletters}{"0A}

%% Math accents
\DeclareMathAccent{\acute} {\mathalpha}{cmoperators}{'023}
\DeclareMathAccent{\bar}   {\mathalpha}{cmoperators}{'026}
\DeclareMathAccent{\breve} {\mathalpha}{cmoperators}{'025}
\DeclareMathAccent{\check} {\mathalpha}{cmoperators}{'024}
\DeclareMathAccent{\ddot}  {\mathalpha}{cmitoperators}{'177}
\DeclareMathAccent{\dot}   {\mathalpha}{cmitoperators}{'137}
\DeclareMathAccent{\grave} {\mathalpha}{cmoperators}{'022}
\DeclareMathAccent{\hat}   {\mathalpha}{cmoperators}{'136}
\DeclareMathAccent{\tild}  {\mathalpha}{cmoperators}{'176}
\DeclareMathAccent{\vec}   {\mathord}  {cmletters}  {"7E}

\DeclareMathSymbol{\barbar}{\mathalpha}{cmoperators}{'026}

\definecolor{titleblue}{RGB}{173,216,230}

%% ---- Page geometry --------------------------------------------------------
\linespread{1.5}
\setlength{\parindent}{0pt}
\setlength{\oddsidemargin}{0pt}
\setlength{\voffset}{-30pt}
\setlength{\textwidth}{450pt}
\setlength{\textheight}{650pt}
\setlength{\headheight}{13.6pt}

%% ---- Section numbering style: "1." / "1.1" --------------------------------
\renewcommand{\thesection}{\@arabic\c@section.}
\renewcommand{\thesubsection}{\thesection\@arabic\c@subsection}

%% ---- Internal storage for title / author / date ---------------------------
\newcommand{\TheTitle}{}
\newcommand{\TheAuthor}{}
\newcommand{\TheDate}{}

\renewcommand{\title}[1]{\renewcommand{\TheTitle}{#1}}
\renewcommand{\author}[1]{\renewcommand{\TheAuthor}{#1}}
\renewcommand{\date}[1]{\renewcommand{\TheDate}{#1}}

%% ---- fancyhdr: running header and footer ----------------------------------
\pagestyle{fancy}
\lhead{}
\chead{}
\cfoot{}
\lfoot{\bfseries\thepage~/ \pageref{LastPage}}
\rfoot{\bfseries\TheAuthor~/ \TheDate}
\rhead{\bfseries\TheTitle}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}

%% Plain style used on the title page (no header rule, but keep footer rule)
\fancypagestyle{plain}{%
  \fancyhf{}%
  \lfoot{\bfseries\thepage~/ \pageref{LastPage}}%
  \rfoot{\bfseries\TheAuthor~/ \TheDate}%
  \renewcommand{\headrulewidth}{0pt}%
  \renewcommand{\footrulewidth}{0.4pt}}

%% ---- Boxed title ----------------------------------------------------------
\renewcommand{\maketitle}{%
  \thispagestyle{plain}%
  \vspace*{-40pt}%
  \begin{tabular}{|c|}\hline
    \rowcolor{titleblue}%
    \parbox{15.46cm}{\vspace*{10pt}\centering\LARGE{\bfseries\TheTitle}\par\vspace{10pt}}%
  \\ \hline
  \end{tabular}%
  \par\vspace{6pt}%
}
