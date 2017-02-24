@echo off
if [%1]==[] goto usage

set _CYGBIN=%1\bin
if not exist "%_CYGBIN%" echo Could not find Cygwin at "%_CYGBIN%" & exit 3

chdir ..
set CHERE_INVOKING=1
%_CYGBIN%\bash --login -c 'sh Make.sh style'
%_CYGBIN%\bash --login -c 'sh Make.sh models'

@echo Done.
goto :eof
:usage
@echo Usage: %0 ^<CygwinPath^>
exit /B 1