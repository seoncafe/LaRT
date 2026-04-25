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
%% \title, \author, \date behave exactly as in the original kiseon_memo.tex:
%%   they store values used in the header/footer and the boxed \maketitle.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{kiseon_memo}[2026/04/18 v1.0 Kiseon memo class]

%% ---- Pass all unknown options to 'article' --------------------------------
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass{article}

%% ---- Required packages ----------------------------------------------------
\RequirePackage[T1]{fontenc}
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
