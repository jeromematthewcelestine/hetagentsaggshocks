function aggs = compute_aggregates(opt, params, dist, kp, x)

output_mesh = x .* opt.z_mesh .* (opt.k_mesh.^params.alpha);

aggs = struct();

aggs.output		= sum(sum(output_mesh .* reshape(dist,opt.n_k,opt.n_z)));
investment_mesh = kp - (1-params.delta).* opt.k_mesh;
aggs.investment = sum(sum(investment_mesh .* reshape(dist,opt.n_k,opt.n_z)));
consumption_mesh = output_mesh - investment_mesh;
aggs.consumption = sum(sum(consumption_mesh .* reshape(dist,opt.n_k,opt.n_z)));
