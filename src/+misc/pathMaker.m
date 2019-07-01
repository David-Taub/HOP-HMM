function ret = pathMaker(params, N, L, prefix, suffix)
    ret = sprintf('%s_m%dbg%dk%do%db%ds%dr%dN%dL%d%s', prefix, params.m, params.backgroundAmount, params.k, ...
                  params.order, params.doGTBound, params.doESharing, params.doResampling, N, L, suffix)
end