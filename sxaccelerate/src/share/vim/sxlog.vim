" Syntax folding for SxAccelerate debug output
"
" 1. copy to ~/.vim/syntax/sxlog.vim
" 2. add to  ~/.vim/syntax/syntax.vim
" au BufRead *.sxlog so $HOME/.vim/syntax/sxlog.vim
"
syn clear

syn region sxBlock  matchgroup=sxBlockStart
                  \ start="^___\([^:]*:[0-9]* [^:]*\): BEGIN$"
                  \ end=": END"
                  \ fold transparent extend


set foldmethod=syntax


set foldtext=SxBlockText()
function! SxBlockText()
  let line = getline(v:foldend)
  let sub = substitute (line, '^___\([^:]*:[0-9]* [^:]*\): END \([0-9.]*\) ms', '\2ms: \1 ', 'g')
  let nNestedBlocks=0
  let iLine = v:foldstart + 1
  while iLine < v:foldend
     let line = getline (iLine)
     if (line =~ ': BEGIN')
        let nNestedBlocks = nNestedBlocks + 1
     endif
     let iLine = iLine + 1
  endwhile
  let sub = "[" . nNestedBlocks . "] " . sub


  let n = v:foldend - v:foldstart + 1
  let info = " " . n . " lines"
  let sub = sub . "                                                                                                                  "
  let num_w = getwinvar( 0, '&number' ) * getwinvar( 0, '&numberwidth' )
  let fold_w = getwinvar( 0, '&foldcolumn' )
  let sub = strpart( sub, 0, winwidth(0) - strlen( info ) - num_w - fold_w - 1 )
  return sub . info
endfunction


hi link sxBlock Statement
