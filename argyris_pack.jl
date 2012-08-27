# ------------------------------------------------------------------------------
# Load the shared object file.
# ------------------------------------------------------------------------------
libap = dlopen("libargyris_pack.so")

__ap_local_functions   = dlsym(libap, :ap_local_functions)
__ap_local_gradients   = dlsym(libap, :ap_local_gradients)
__ap_local_hessians    = dlsym(libap, :ap_local_hessians)
__ap_global_maps       = dlsym(libap, :ap_global_maps)
__ap_global_functions  = dlsym(libap, :ap_global_functions)
__ap_global_gradients  = dlsym(libap, :ap_global_gradients)
__ap_global_hessians   = dlsym(libap, :ap_global_hessians)
__ap_matrix_mass       = dlsym(libap, :ap_matrix_mass)
__ap_matrix_stiffness  = dlsym(libap, :ap_matrix_stiffness)
__ap_matrix_biharmonic = dlsym(libap, :ap_matrix_biharmonic)

# ------------------------------------------------------------------------------
# Julia interfaces to the .so file.
# ------------------------------------------------------------------------------
function ap_local_functions{T}(x::Vector{T}, y::Vector{T})
# evaluate the 21 Argyris basis functions at quadrature points specified by
# 'points'.
    check_size(x,y)

    function_values = zeros(21,size(x)[1])
    ccall(__ap_local_functions, Void,
          (Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}),
          x, y, length(x), function_values)

    return function_values
end

function ap_local_gradients{T}(x::Vector{T}, y::Vector{T})
# evaluate the derivatives of the 21 Argyris basis functions at quadrature
# quadrature points specified by 'points'.
    check_size(x,y)

    dx = zeros(21,size(x)[1])
    dy = zeros(21,size(x)[1])
    ccall(__ap_local_gradients, Void,
          (Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}, Ptr{Float64}),
          x, y, size(x)[1], dx, dy)

    return dx, dy
end

function ap_local_hessians{T}(x::Vector{T}, y::Vector{T})
# evaluate the second derivatives of the 21 Argyris basis functions at quadrature
# quadrature points specified by 'points'.
    check_size(x,y)

    dxx = zeros(21,size(x)[1])
    dxy = zeros(21,size(x)[1])
    dyy = zeros(21,size(x)[1])
    ccall(__ap_local_hessians, Void,
          (Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}),
          x, y, length(x), dxx, dxy, dyy)

    return dxx, dxy, dyy
end

function ap_global_maps{T}(x::Vector{T}, y::Vector{T})
# Use the change-of-basis matrix C to convert reference values to global values.
    check_size(x,y)
    if length(x) != 3
        error("There should be three corners to a triangle.")
    end

    C = zeros(21,21)
    B = zeros(2,2)
    b = zeros(2)
    Th = zeros(3,3)
    ccall(__ap_global_maps, Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}),
          x, y, C, B, b, Th)

    return C, B, b, Th
end

function ap_global_functions(C, ref_values)
# Use the change-of-basis matrix C to convert reference values to global values.
    check_transformations(C)
    check_reference_values(ref_values)

    function_values = zeros(21,size(ref_values)[2])
    ccall(__ap_global_functions, Void, (Ptr{Float64}, Ptr{Float64}, Int64,
          Ptr{Float64}), C, ref_values, size(ref_values)[2],
          function_values)

    return function_values
end

function ap_global_gradients(C, B, ref_dx, ref_dy)
# evaluate the derivatives of the 21 Argyris basis functions given values of
# the derivatives on the reference triangle.
    check_transformations(C,B)
    check_reference_values(ref_dx, ref_dy)

    dx = zeros(size(ref_dx))
    dy = zeros(size(ref_dx))
    ccall(__ap_global_gradients, Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64,
           Ptr{Float64}, Ptr{Float64}),
          C, B, ref_dx, ref_dy, size(ref_dx)[2], dx, dy)

    return dx, dy
end

function ap_global_hessians(C, Th, ref_dxx, ref_dxy, ref_dyy)
# Evaluate second derivatives of the basis functions given reference values.
    check_transformations(C)
    if size(Th) != (3,3)
        error("Incorrect dimensions in the 3x3 'Theta' matrix.")
    end
    check_reference_values(ref_dxx, ref_dxy, ref_dyy)

    dxx = zeros(size(ref_dxx))
    dxy = zeros(size(ref_dxx))
    dyy = zeros(size(ref_dxx))
    ccall(__ap_global_hessians, Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          C, Th, ref_dxx, ref_dxy, ref_dyy, size(ref_dxx)[2], dxx, dxy, dyy)

    return dxx, dxy, dyy
end

function ap_matrix_mass(C, B, ref_functions, weights)
# Evaluate the local mass matrix given reference values.
    check_transformations(C,B)
    check_reference_values(ref_functions)
    check_weights(ref_functions, weights)

    mass = zeros(21,21)
    ccall(__ap_matrix_mass, Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64,
           Ptr{Float64}),
          C, B, ref_functions, weights, size(weights)[1], mass)

    return mass
end

function ap_matrix_stiffness(C, B, ref_dx, ref_dy, weights)
# Evaluate the local stiffness matrix given reference values.
    check_transformations(C,B)
    check_reference_values(ref_dx, ref_dy)
    check_weights(ref_dx, weights)

    stiffness = zeros(21,21)
    ccall(__ap_matrix_stiffness, Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Int64, Ptr{Float64}),
          C, B, ref_dx, ref_dy, weights, size(weights)[1], stiffness)

    return stiffness
end

function ap_matrix_biharmonic(C, B, Th, ref_dxx, ref_dxy, ref_dyy, weights)
# Evaluate the local biharmonic matrix given reference values.
    check_transformations(C,B,Th)
    check_reference_values(ref_dxx, ref_dxy, ref_dyy)
    check_weights(ref_dxx, weights)

    biharmonic = zeros(21,21)
    ccall(__ap_matrix_biharmonic, Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}),
          C, B, Th, ref_dxx, ref_dxy, ref_dyy, weights, size(weights)[1],
          biharmonic)

    return biharmonic
end

# ------------------------------------------------------------------------------
# These functions pertain to validating input.
# ------------------------------------------------------------------------------
function check_size(arrays...)
# raise an exception if the arrays are not equal in size.
    if length(arrays) == 1
        return
    else
        for i=2:length(arrays)
            if  size(arrays[i-1]) != size(arrays[i])
                error("Mismatch in array sizes.")
            end
        end
    end
end

function check_transformations(transformations...)
# check the sizes of the matrix transformations C, B, and Th.
    if size(transformations[1]) != (21,21)
        error("Incorrect dimensions in the 21x21 'C' matrix.")
    end
    if length(transformations) > 1
        if size(transformations[2]) != (2,2)
            error("Incorrect dimensions in the 2x2 'B' matrix.")
        end
    end
    if length(transformations) > 2
        if size(transformations[3]) != (3,3)
            error("Incorrect dimensions in the 3x3 'Theta' matrix.")
        end
    end
end

function check_reference_values(values...)
# Raise an exception if the input arrays do not all have 21 rows and the same
# number of columns.
    if size(values[1])[1] != 21
        error("There should be 21 rows corresponding to 21 basis functions.")
    end
    if length(values) == 1
        return
    else
        for i=2:length(values)
            if  size(values[i-1]) != size(values[i])
                error("Mismatch in array sizes.")
            end
        end
    end
end

function check_weights{T}(values, weights::Vector{T})
    if size(weights)[1] != size(values)[2]
        error("Mismatch in number of weights and number of quadrature points.")
    end
end