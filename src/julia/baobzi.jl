module baobzi

mutable struct baobzi_t
end

mutable struct baobzi_input_t
    func::Ptr{Cvoid}
    data::Ptr{Cvoid}
    dim::Cint
    order::Cint
    tol::Cdouble
    minimum_leaf_fraction::Cdouble
    split_multi_eval::Cint
end

function init(fin, dim, order, center, half_length, tol, minimum_leaf_fraction=0.0, split_multi_eval=1)
    fanon = (x, p) -> fin(x)
    fbind = @cfunction($fanon, Cdouble, (Ptr{Cdouble}, Ptr{Cvoid},))
    input = baobzi_input_t(fbind.ptr, C_NULL, dim, order, tol, minimum_leaf_fraction, split_multi_eval)

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

function eval_multi(func::Ptr{baobzi_t}, x)
    n_points = Int32(length(x) // 2)
    println(n_points)
    res = Array{Float64}(undef, n_points)
    ccall(
        (:baobzi_eval_multi, :libbaobzi),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Cint,),
        func, x, res, n_points
    )
    return res
end

function stats(func::Ptr{baobzi_t})
    ccall(
        (:baobzi_stats, :libbaobzi),
        Ptr{Cvoid},
        (Ptr{Cvoid},),
        func
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
