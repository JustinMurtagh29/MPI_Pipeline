if ispc
    % Windows7/Visual Studio 2010
    mex -ljson -g LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:MSVSRT.lib /NODEFAULTLIB:LIBCMT.lib" -IC:\Users\panton\code\ -LC:\Users\panton\code\json\Release\ fromjson.c
    mex -ljson -g LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:MSVSRT.lib /NODEFAULTLIB:LIBCMT.lib" -IC:\Users\panton\code\ -LC:\Users\panton\code\json\Release\ tojson.c 
else
    % Linux/Ubuntu/GCC
    mex -ljson-c -I/zdata/manuel/code/active/matlab-json -L/zdata/manuel/usr/local/lib -g fromjson.c
    mex -ljson-c -I/zdata/manuel/code/active/matlab-json -L/zdata/manuel/usr/local/lib -lm -g tojson.c
end

mex setjsonfield.c
