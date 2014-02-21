
@ECHO ON
@ECHO %CD%

IF NOT EXIST "C:\Program Files\rpn" MKDIR "C:\Program Files\rpn"
IF NOT EXIST "C:\Program Files\rpn\include" MKDIR "C:\Program Files\rpn\include"
IF NOT EXIST "C:\Program Files\rpn\bin" MKDIR "C:\Program Files\rpn\bin"
IF NOT EXIST "C:\Program Files\rpn\lib" MKDIR "C:\Program Files\rpn\lib"

COPY ..\src\*.h "C:\Program Files\rpn\include"
COPY Debug\*.exe "C:\Program Files\rpn\bin"
COPY Debug\*.lib "C:\Program Files\rpn\lib"
COPY ..\rpn.ico "C:\Program Files\rpn\bin"

PAUSE
