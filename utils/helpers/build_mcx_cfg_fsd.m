function mcx_cfg_full = build_mcx_cfg_fsd(mcx_cfg_base, CaseNow)
% BUILD_MCX_CFG_FSD  Build an MCX configuration for a full source-detector run.
%
%  mcx_cfg_full = build_mcx_cfg_fsd(mcx_cfg_base, CaseNow)
%
%  Branches off the shared base configuration and sets the source type and
%  shape parameters appropriate for the geometry described in CaseNow.
%  All other MCX fields (photon count, optical properties, GPU, etc.) are
%  inherited unchanged from mcx_cfg_base.
%
%  This function encapsulates all the MCX source-type logic so that the
%  main loop remains readable. To add a new geometry type, only this
%  function needs to be extended.
%
%  INPUTS
%  ------
%  mcx_cfg_base  - struct: shared MCX configuration (see Section 2 of the
%                  main demo script). Contains nphoton, vol, prop, srcpos,
%                  unitinmm, boundary conditions, GPU settings, etc.
%                  See mcxlab.m for a description of all fields.
%
%  CaseNow       - struct: geometry descriptor for one scenario, as
%                  returned by build_sia_cases.m. Relevant fields:
%                    .type           'uniform_shapes' | 'square_square'
%                                    | 'gaussian_overlap'
%                    .source_shape   'disk' | 'annulus' | 'square'
%                    .source_rin     inner source radius [mm] (disk/annulus)
%                    .source_rout    outer source radius [mm] (disk/annulus)
%                    .source_side    side length [mm] (square only)
%                    .w_mm           Gaussian waist [mm] (gaussian_overlap)
%
%  OUTPUT
%  ------
%  mcx_cfg_full  - struct: MCX configuration ready to pass to mcxlab().
%                  Differs from mcx_cfg_base only in the source fields:
%                    .srctype    source type string ('disk','planar','gaussian')
%                    .srcparam1  source shape parameter 1 (type-dependent)
%                    .srcparam2  source shape parameter 2 (type-dependent)
%                    .srcpos     (adjusted for planar sources)
%                    .srcdir     (adjusted for Gaussian sources)
%
%  NOTES ON UNITS
%  --------------
%  MCX expects all source size parameters in VOXEL units, not in mm.
%  This function divides all physical dimensions [mm] by unitinmm
%  [mm/voxel] before storing them in the config struct.
%
%  MCX SOURCE TYPE REFERENCE (abbreviated; see mcxlab.m for full details)
%  -----------------------------------------------------------------------
%   'disk'     - circular disk source. srcparam1 = [r_outer, r_inner, 0, 0]
%                in voxels (r_inner = 0 for a full disk, > 0 for an annulus).
%   'planar'   - rectangular (planar) source. srcpos is the corner of the
%                rectangle; srcparam1 = [dx 0 0 0] and srcparam2 = [0 dy 0 0]
%                define the two spanning edge vectors in voxels.
%   'gaussian' - Gaussian beam. srcparam1 = [w 0 0 0] where w is the 1/e^2
%                waist radius in voxels. srcdir = [nx ny nz w_focus] where
%                the 4th element encodes the focus distance (0 = collimated).

% Start from the shared base configuration so that only the source-specific
% fields need to be modified below.
mcx_cfg_full = mcx_cfg_base;

% Remove any leftover source shape parameters from a previous iteration to
% prevent stale values from accidentally influencing the new simulation.
mcx_cfg_full = rmfield_if_exists(mcx_cfg_full, 'srcparam1');
mcx_cfg_full = rmfield_if_exists(mcx_cfg_full, 'srcparam2');

switch CaseNow.type

    % ------------------------------------------------------------------
    % Uniform disk or annular source
    % ------------------------------------------------------------------
    case 'uniform_shapes'
        % MCX 'disk' source type handles both full disks and annuli via
        % srcparam1 = [r_outer, r_inner, 0, 0] in voxel units.
        % r_inner = 0 gives a full disk; r_inner > 0 gives a ring/annulus.
        mcx_cfg_full.srctype = 'disk';

        if strcmp(CaseNow.source_shape, 'disk')
            % Full disk: outer radius only.
            mcx_cfg_full.srcparam1 = [CaseNow.source_rout, 0, 0, 0] / mcx_cfg_full.unitinmm;

        elseif strcmp(CaseNow.source_shape, 'annulus')
            % Annulus: outer and inner radii.
            mcx_cfg_full.srcparam1 = [CaseNow.source_rout, CaseNow.source_rin, 0, 0] / mcx_cfg_full.unitinmm;

        else
            error('build_mcx_cfg_fsd: unsupported source_shape ''%s'' for type ''uniform_shapes''.', ...
                  CaseNow.source_shape);
        end

        % srcparam2 is not used for the disk source type; MCX ignores it,
        % but we set it to zero for clarity.
        mcx_cfg_full.srcparam2 = [0 0 0 0];

    % ------------------------------------------------------------------
    % Square (planar) source
    % ------------------------------------------------------------------
    case 'square_square'
        % MCX 'planar' source type launches photons uniformly from a
        % parallelogram defined by a corner point and two edge vectors.
        % For an axis-aligned square of side `side_vox` centred on srcpos,
        % the corner is at srcpos - [side_vox/2, side_vox/2, 0] and the
        % edge vectors are [side_vox 0 0 0] and [0 side_vox 0 0].
        mcx_cfg_full.srctype = 'planar';

        % Convert side length from mm to voxels.
        side_vox = CaseNow.source_side / mcx_cfg_full.unitinmm;  % [voxels]

        % Shift srcpos to the lower-left corner of the square, keeping
        % the square centred on the original source centre.
        corner_vox = [double(mcx_cfg_base.srcpos(1)) - side_vox/2, ...
                      double(mcx_cfg_base.srcpos(2)) - side_vox/2, ...
                      double(mcx_cfg_base.srcpos(3))];

        mcx_cfg_full.srcpos    = corner_vox;
        % Edge vectors of the source square: x-direction and y-direction.
        mcx_cfg_full.srcparam1 = [side_vox 0 0 0];  % x edge vector [voxels]
        mcx_cfg_full.srcparam2 = [0 side_vox 0 0];  % y edge vector [voxels]
        mcx_cfg_full.srcdir    = [0 0 1];            % launch direction: +z (into medium)

    % ------------------------------------------------------------------
    % Gaussian beam source
    % ------------------------------------------------------------------
    case 'gaussian_overlap'
        % MCX 'gaussian' source type launches photons with a Gaussian
        % spatial intensity profile. srcparam1(1) is the 1/e^2 beam waist
        % radius in voxel units. The 4th element of srcdir encodes the
        % focal distance; 0 means collimated (parallel beam at the surface).
        mcx_cfg_full.srctype   = 'gaussian';
        mcx_cfg_full.srcdir    = [0 0 1 0];  % [nx ny nz focus_distance_vox]
        mcx_cfg_full.srcparam1 = [CaseNow.w_mm / mcx_cfg_full.unitinmm, 0, 0, 0];  % [waist_vox 0 0 0]
        mcx_cfg_full.srcparam2 = [0 0 0 0];  % unused for Gaussian source

    otherwise
        error('build_mcx_cfg_fsd: unknown case type ''%s''.', CaseNow.type);
end

end % function build_mcx_cfg_fsd

% -------------------------------------------------------------------------
% Local helper
% -------------------------------------------------------------------------
function S = rmfield_if_exists(S, fieldname)
    % Remove a field from struct S only if it is present.
    % Avoids errors when the base config does not have the field.
    if isfield(S, fieldname)
        S = rmfield(S, fieldname);
    end
end
