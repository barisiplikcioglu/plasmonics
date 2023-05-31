function [eff,sp] = spp_efficiency_lalanne(lambda,w,em,n1,n2,source,angle)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% This function computes the complex SPP generation strength/amplitude and
% reflection/transmission coefficient at a subwavelength, symmetric 2D metallic 
% slit excited either by a plane wave or the fundamental slit mode. The code
% uses the approximate mode-matching method by P. Lalanne et al., which is 
% detailed in the following publications:
% 
% * P. Lalanne, J. P. Hugonin, and J. C. Rodier, "Theory of surface plasmon 
%   generation at nanoslit apertures". Phys. Rev. Lett. 95. p. 263902. (2005)
% * P. Lalanne, J. P. Hugonin, and J. C. Rodier, "Approximate model for 
%   surface-plasmon generation at slit apertures," J. Opt. Soc. Am. A 23,
%   p. 1608-1615 (2006)
% 
% INPUTS:
% * lambda: Free-space wavelength of the light (same units as slit width 'w')
% * w: Slit width (same units as free-space wavelength 'lambda')
% * em: Complex relative permittivity of the metal
% * n1: Refractive index of the dielectric inside the slit (lossless)
% * n2: Refractive index of the dielectric outside the slit (lossless)
% * source: Excitation source, 'pw' for plane wave or 'slit' for the slit mode
% * angle: Incidence angle with respect to slit normal for the plane wave 
%          source in degrees, if 'pw' option is also chosen. This parameter
%          may be left blank for slit excitation.
%
% OUTPUTS:
% * eff: Complex SPP amplitude (defined as SPP generation strength) at slit 
%        for the chosen excitation source. SPP excitation efficiency can be
%        computed from the modulus squared of this value (2*|eff|^2).
% * sp: Scattering parameter for the fundamental slit mode. Depending on
%       the excitation source, it corresponds to the following parameters:
%       * For slit excitation ('slit'), 'sp' corresponds to the reflection
%         coefficient for the fundamental mode that impinges on the slit.
%       * For plane wave excitation ('pw'), 'sp' corresponds to the 
%         transmission coefficient of the fundamental slit mode excited by
%         the impinging plane wave.
% 
% EXAMPLE SYNTAX:
% spp_efficiency_lalanne(800e-9,80e-9,-26.27+1.85i,1,1,'pw',15)
% spp_efficiency_lalanne(800e-9,80e-9,-26.27+1.85i,1,1,'slit')
%
% LIMITATIONS:
% The approximate model is best suited for longer wavelengths (at which the
% metal approximately behaves like PEC), as well as for subwavelength slits 
% and non-oblique plane wave incidence (angle<=30), for which single mode 
% approximation holds. The validity and limitations of the model is outlined
% in section 4C of the 2006 paper.
%
% NOTES:
% * This formulation uses exp(-i*omega*t) convention.
% * In Lalanne et al.'s original formulation, the reflection coefficients are
%   defined with respect to the transverse magnetic field of the fundamental 
%   mode, rather than the transverse electric field. Thus for all complex results, 
%   phase data must be interpreted accordingly.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if nargin > 7
    error ('Too many arguments.')
elseif ((nargin==6) && strcmp(source,'pw'))
    error ('No incidence angle is defined for the plane wave source.')
elseif nargin < 6
    error('Too few arguments.')
end

if nargin==6
    angle=0;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Note: The following function computes the square root according to 
%       Sommerfeld radiation condition.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
sqrtbranch=@(x) sqrt(x).*(-1).^(imag(x)<0);

sincx=@(x) sin(x)./x;
weff=w./(lambda./n2);
i0fun=@(u) sincx(pi*weff*u).^2./sqrtbranch(1-u.^2);
i0=quadgk(i0fun,-Inf,-1)+quadgk(i0fun,-1,1)+quadgk(i0fun,1,Inf);
i1fun=@(u) sincx(pi*weff.*u).*exp(-1i*weff*pi.*u)./...
           sqrtbranch(1-u.^2)./(sqrtbranch(1-u.^2)+sqrtbranch(n2.^2./(em+n2.^2)));
i1=quadgk(i1fun,-Inf,-1)+quadgk(i1fun,-1,1)+quadgk(i1fun,1,Inf);

ref=(n2*weff*i0-n1)./(n2*weff*i0+n1);
alpha=-1i*sqrtbranch(4./pi*n2.^2./n1.*sqrt(abs(em))./(-em-n2.^2).*weff).*...
            i1./(1+(n2./n1).*weff.*i0);  

switch source
    case 'slit'
        sp=ref;
        eff=alpha;                                                         
    case 'pw'
        N0=w./n1;
        NP=w.*cosd(angle)./n2;
        sp=sqrt(N0./NP).*2.*sinc(weff*sind(angle))./((n2./n1).*weff.*i0+1);
        eff=-sqrt(N0./NP).*sinc(weff.*sind(angle)).*alpha;
    otherwise
        error('Please enter a valid excitation source.')
end

end
