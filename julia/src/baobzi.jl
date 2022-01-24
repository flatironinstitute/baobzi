module baobzi

mutable struct baobzi_t
end

mutable struct baobzi_input_t
    func::Ptr{Cvoid}
    data::Ptr{Cvoid}
    dim::Cint
    order::Cint
    tol::Cdouble
end

function init(fin, dim, order, center, half_length, tol)
    fanon = (x, p) -> fin(x)
    fbind = @cfunction($fanon, Cdouble, (Ptr{Cdouble}, Ptr{Cvoid},))
    input = baobzi_input_t(fbind.ptr, C_NULL, dim, order, tol)

    output_ptr = ccall(
        (:baobzi_init, :libbaobzi),
        Ptr{baobzi_t},
        (Ref{baobzi_input_t}, Ptr{Cdouble}, Ptr{Cdouble},),
        input, center, half_length
    )

    return output_ptr
end

function eval(func::Ptr{baobzi_t}, x)
    return ccall(
        (:baobzi_eval, :libbaobzi),
        Cdouble,
        (Ptr{Cvoid}, Ptr{Cdouble},),
        func, x
    )
end

function free(func::Ptr{baobzi_t})
    ccall(
        (:baobzi_free, :libbaobzi),
        Ptr{Cvoid},
        (Ptr{Cvoid},),
        func
    )
end

function save(func::Ptr{baobzi_t}, filename::String)
    ccall(
        (:baobzi_save, :libbaobzi),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Cstring),
        func, filename
    )
end

function restore(filename::String)
    output_ptr = ccall(
        (:baobzi_restore, :libbaobzi),
        Ptr{baobzi_t},
        (Cstring,),
        filename
    )

    return output_ptr
end

end
