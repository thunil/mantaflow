REM parameters: 		                     openB bWidth  sizeFactor  
..\build\Release\manta.exe plume_2d.py        yY     1      1 1 1  
    

REM todo zZ open ..\build\Release\manta.exe plume_2d.py       xXyY     1    1    1 1 1         0   

REM finish
REM ..\build\Release\manta.exe plume_2d.py     xXyY     1    1    4 4 1         0   
REM ..\build\Release\manta.exe plume_2d.py     xXyY     1    8    4 4 1         0   
REM ..\build\Release\manta.exe plume_2d.py     xXyY     0    4    4 4 1         0   
REM ..\build\Release\manta.exe plume_2d.py     xXyY     0    1    4 4 1         0   
REM ..\build\Release\manta.exe plume_2d.py     xXyY     0    8    4 4 1         0   

REM now all fluid for some xX
REM now only extra from fluid cell -> check empty, half fluid, full fluid

REM now levelset, flip evtl high bulkvel -> unterschied testen mit und ohne eigener extrapolation
REM braucht eig keine Analyse der verschiedenen Parametern, auch kein inflow da nie inflow (oder??)
REM keine vortices bei flip, levelset -> schwer richtige extrapolation zu testen
