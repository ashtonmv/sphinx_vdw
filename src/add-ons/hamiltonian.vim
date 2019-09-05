" Vim syntax file
" Language:	multiband k*p-Hamiltonians
" Maintainer:	Oliver Marquardt; o.marquardt@yahoo.ie
" Last Change:	26.1.2012
" Version:	1
" How to get the syntax highlighting running:
"
" put the following lines in your .vimrc which should be directly in the home folder.
" If .vimrc is not there, create it.
"
" 	so $HOME/.vim/syntax/syntax.vim
" 	sy on
"
" in ./vim/syntax/syntax.vim put the following:
" 	augroup filetypedetect
" 	au! BufRead,BufNewFile *.ham	setfiletype hamiltonian
" 	augroup END
"
" copy this file (hamiltonian.vim) into the ~/.vim/syntax/ folder.
" If you open any *.ham-file, the syntax highlighting should be active now.

syn clear

syn keyword reserved kx ky kz k kx2 ky2 kz2 k2 kxy kxz kyz
syn keyword reserved eX eY eZ eXX eYY eZZ eXY eXZ eYZ e
syn keyword reserved Vp Vext Hamiltonian
syn keyword maths Sqrt i
syn match  sxNumber "[-+]\=\(\<\d[[:digit:]_]*L\=\>\|0[xX]\x[[:xdigit:]_]*\>\)"

" comments follow mathematica style
syn region mathComment start=+/\*+ end=+\*/+ skipempty
syn region sqrtBracket start=+{+ end=+}+ skipempty

hi link reserved Type
hi link maths Include
hi link sxNumber Include
hi link sqrtBracket Include
hi link mathComment       Comment

let b:current_syntax = "hamiltonian"
