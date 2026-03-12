function a = input_type(A)

  % case 1: Function Handle
  if isa(A, "function_handle")
    a = 1;
    return;
  end

  % case 2: Struct with field 'E'
  if isstruct(A) && (isfield(A, "E") || isfield(A, "Ep"))
    a = 2;
    return;
  end

  % case 3: anything else (Q)
  a = 3;

end
