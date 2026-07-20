using RetrievalToolbox

@info "Quick check to see if XRTM communicates with RetrievalToolbox..."

max_coef = 10
n_quad = 2
n_stokes = 1
n_derivs = 0
n_layers = 10
n_theta_0s = 1
n_kernel_quad = 32
n_out_levels = 1
n_out_thetas = 1

@info "Creating XRTM object .. "

xrtm_obj = RetrievalToolbox.XRTM.create(
    ["source_thermal", "output_at_levels"],
    ["eig_bvp"],
    max_coef,
    n_quad,
    n_stokes,
    n_derivs,
    n_layers,
    n_theta_0s,
    n_kernel_quad,
    ["lambertian"],
    n_out_levels,
    n_out_thetas
)

@info "XRTM object is $(xrtm_obj)"

@info "Check max. coefs"
@assert RetrievalToolbox.XRTM.get_max_coef(xrtm_obj) == max_coef

@info "Check N quadratures"
@assert RetrievalToolbox.XRTM.get_n_quad(xrtm_obj) == n_quad

@info "Check N Stokes components"
@assert RetrievalToolbox.XRTM.get_n_stokes(xrtm_obj) == n_stokes

@info "Check number of derivatives"
@assert RetrievalToolbox.XRTM.get_n_derivs(xrtm_obj) == n_derivs

@info "Check number of layers"
@assert RetrievalToolbox.XRTM.get_n_layers(xrtm_obj) == n_layers

@info "Check number of surface kernel quadratures"
@assert RetrievalToolbox.XRTM.get_n_kernel_quad(xrtm_obj) == n_kernel_quad

@info "Check number of surface kernels"
@assert RetrievalToolbox.XRTM.get_n_kernels(xrtm_obj) == 1

@info "Check number of outgoing levels"
@assert RetrievalToolbox.XRTM.get_n_out_levels(xrtm_obj) == n_out_levels

@info "Check number of outgoing viewing zenith angles"
@assert RetrievalToolbox.XRTM.get_n_out_thetas(xrtm_obj) == n_out_thetas

@info "Destroying XRTM object"
RetrievalToolbox.XRTM.destroy(xrtm_obj)


@info "All done!"