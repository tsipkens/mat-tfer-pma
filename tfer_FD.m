
% TFER_FD   Evaluates the transfer function of a PMA using finite differences.
% Author:   Timothy Sipkens, 2018-12-27
% 
% Inputs:
%   sp          Structure defining various setpoint parameters (e.g. m_star, V)
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length)
%
% Outputs:
%   Lambda      Transfer function
%   sp          Device properties (e.g. classifier length)
%   n           Struct containing information about the particle
%               distribution at different axial position (used for plotting)
%
% Notes:
% 	This code is adapted from Buckley et al. (2017) and the Olfert
% 	laboratory. 
%=========================================================================%

function [Lambda,n] = tfer_FD(sp,m,d,z,prop)

[tau, C0, D] = parse_inputs(sp,m,d,z,prop);
        % parse inputs for common parameters


%-- Discretize the space -------------------------------------------------%
nr = 200;
dr = (prop.r2 - prop.r1) / nr;
r_vec = (prop.r1 + dr/2):dr:(prop.r2 - dr/2);

nz = 200;
dz = prop.L /( nz-1);
if nargout >= 2; n.z_vec = 0:dz:prop.L; n.r_vec = r_vec; end % vector of z positions, used for plotting only


%-- Parameters related to CPMA geometry ----------------------------------%
D0 = D .* prop.L / (prop.del ^ 2 * prop.v_bar);
	% dimensionless diffusion coefficient
    

%-- Evaluate relevant radial quantities ----------------------------------%
v_theta = [sp.alpha]' .* r_vec + [sp.beta]' ./ r_vec; % azimuthal velocity distribution
F_e = C0 ./ r_vec; % electrostatic force
v_z = 3/2 .* prop.v_bar .* (1 - ((r_vec - prop.rc) ./ (prop.del)) .^ 2);
    % axial velocity distribution (parabolic)
% v_z = prop.v_bar.*ones(size(r_vec)); % axial velocity distribution (plug)


Lambda = zeros(length(sp), length(m));% initialize the transfer function variable
for jj=1:length(sp)
    %-- Speed computation using resolution to limit computation --------------%
    ind = 1:length(m);
    if isfield(sp(jj), 'Rm') % if resolution is specified, use to reduce necessary computation
        cond0 = or(m > (z .* sp(jj).m_star + 2 .* sp(jj).m_max),...
            m < (z .* sp(jj).m_star - 2 .* sp(jj).m_max));
                % NOTE: conditions limits consideration of those regions where particles
                % do not escape, speeding computation.
        ind = ind(~cond0);
    end


    %-- Loop over mass, but not m_star -----------------------------------%
    for ii=ind % loop over mass (not m_star)

        F_c = m(ii) .* v_theta(jj,:) .^ 2 ./ r_vec; % centriputal force
        c_r = tau(ii) / m(ii) .* (F_c - F_e(jj,:)); % particle velocity across streamlines
        drcr_dr = tau(ii) .* (2 .* sp(jj).alpha .^ 2 .* r_vec - ...
            2 .* sp(jj).beta .^ 2 ./ (r_vec .^ 3));
        
        
        %-- Initialize variables -----------------------------------%
        n_vec = ones(size(r_vec)); % number concentration at current axial position
        n_vec0 = n_vec; % initial number concentration (used for  func. eval.)
        if sp(jj).m_star >= 2 % initilize variables used for visualizing number concentrations
            n_mat = zeros(nz, length(r_vec));
            n_mat(1,:) = n_vec;
        end


        %-- Get coefficients -------------------------------------%
        zet = v_z ./ dz;
        gam = D(ii) / (2 * dr ^ 2);
        kap = D(ii) ./ (4 .* r_vec .* dr);
        eta = 1 ./ (2 .* r_vec) .* drcr_dr;
        phi = c_r ./ (4 * dr);

        %== Crank-Nicholson ========================================%
        [a, b, c, RHS] = gen_eqs(zet, gam, kap, eta, phi);

        %-- Primary loop for finite difference ---------------------%
        for kk = 2:nz
            n_vec = tridiag([0,a], b, c, RHS(n_vec));
                % solve system using Thomas algorithm

            if D0(ii)<1e-3 % for low diffusion, stabalize by smoothing oscillations
                n_vec = [0,n_vec,0];
                n_vec = (n_vec(1:end-2) + 1e2 .* n_vec(2:end-1) + ...
                    n_vec(3:end)) ./ ...
                    (1e2 + 2);
            end

            if nargout>=2; n_mat(kk,:) = n_vec; end
                % record particle distribution at this axial position
        end

        %-- Format output ------------------------------------------%
        if nargout>=2 % for visualizing number concentrations
            n.n_mat{jj,ii} = max(n_mat,0);
        end

        Lambda(jj,ii) = sum(n_vec .* v_z) ./ sum(n_vec0 .* v_z);
            % evaluate transfer fucntion

        if Lambda(jj,ii) < 0.0005; Lambda(jj,ii) = 0; end
            % truncate small  func. values

    end
end

end


%== GEN_EQS ==============================================================%
%   Generate matrix equations for Crank-Nicholson scheme
%   Author:   Timothy Sipkens, 2018-12-27
function [a, b, c, RHS] = gen_eqs(zet, gam, kap, eta, phi)

ind_mid = 2:(length(zet)-1);


%-- Righ-hand side of eq. to be solved --------------------%
RHS1 = @(n) (zet(1) - 2 * gam - eta(1)).*n(1)+...
    (gam+kap(1) - phi(1)) .* n(2);
RHS2 = @(n) (zet(ind_mid) - 2 * gam-eta(ind_mid)) .* n(2:(end-1)) + ...
    (gam - kap(ind_mid) + phi(ind_mid)) .* n(1:(end-2)) + ...
    (gam + kap(ind_mid) - phi(ind_mid)) .* n(3:end);
RHS3 = @(n) (zet(end) - 2 * gam - eta(end)) .* n(end)+...
    (gam - kap(end) + phi(end)) .* n(end-1);
RHS = @(n) [RHS1(n), RHS2(n), RHS3(n)];


%-- Form tridiagonal matrix ------------------------------%
b1 = zet(1) + 2 .* gam + eta(1); % center
c1 = -gam - kap(1) + phi(1); % +1

a2 = (-gam + kap(ind_mid) - phi(ind_mid)); % -1
b2 = zet(ind_mid) + 2 * gam + eta(ind_mid); % center
c2 = -gam - kap(ind_mid) + phi(ind_mid); % +1

a3 = -gam + kap(end) - phi(end); % -1
b3 = zet(end) + 2 * gam+eta(end); % center

a = [a2, a3];
b = [b1, b2, b3];
c = [c1, c2];

end


