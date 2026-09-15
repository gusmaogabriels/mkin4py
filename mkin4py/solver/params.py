"""Convergence controls for the original Newton/RK4 procedure."""
default = {
    'h': 1., 'hfun': 0.995, 'delta_min': 1e-30, 'criteria': 1e-8,
    'inner_criteria': 1e-9, 'convtol': 100, 'convtolH': 20,
    'inner_convtol': 300, 'max_restarts': 100, 'max_time': 60.,
    'linear_solver': 'qmr',
}
convergence_params = default.copy()
