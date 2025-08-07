function PSF = generate_PSF(ap_field, psf_pupil, zernike_pupil, ...
    field_size, lamda, diff_limited, defocus)
%UNTITLED2 Generates PSF for a given field size, pupil size, and test
%wavelength
%   ap_field is the field size, in pixels; typically 512;
%   psf_pupil is the pupil size for testing;
%   zernike_pupil is the pupil size for the wavefront testing (usually == psf_pupil)
%   field_size is the field size, in ARCMIN; (60 arcmin per degree, so 72 for a 1.2deg field)
%   lamda is the stimulus wavelength;
%   diff_limited: if diff_limited == 1, set all Zernike coefficients to
%   zero; else, input a given .zer from HSWS

global PARAMS;	%set as global so that PARAMS(6) can be changed in a subroutine

PARAMS(1) = ap_field; 	% size of pupil aperture field in pixels (this defines the resolution of the calculation)
PARAMS(5) = lamda;	% imaging wavelength in microns
PARAMS(6) = 0;		% number of pixels over which PSF is calculated (do not adjust here set in Zernikephase subroutine)
PARAMS(7) = 20; % increase to enhance the display of the wavefront (doesn't affect calculation)

if diff_limited == 1;
    % Zernike coeffs
    % Preset all coefficients to 0 and edit each value as desired
    
    % tilt
    c(1)=0; 	c(2)=0.0;
    %defocus and astigmatism
    c(3)=0; 	c(4)=0;	c(5)=0;
    %coma like
    c(6)=0.0;	c(7)=0.00;	c(8)=0;	c(9)=0.0;
    %spherical aberration like
    c(10)=0;	c(11)=0.0;	c(12)=00;	c(13)=0;	c(14)=0;
    %higher order (each row is a new radial order)
    c(15)=0;	c(16)=0;	c(17)=0;	c(18)=0;	c(19)=0;	c(20)=0;
    c(21)=0;	c(22)=0;	c(23)=0;	c(24)=0;	c(25)=0;	c(26)=0;	c(27)=0;
    c(28)=0;	c(29)=0;	c(30)=0;	c(31)=0;	c(32)=0;	c(33)=0;	c(34)=0;	c(35)=0;
    c(36)=0;	c(37)=0;	c(38)=0;	c(39)=0;	c(40)=0;	c(41)=0;	c(42)=0;	c(43)=0;	c(44)=0;
    c(45)=0;	c(46)=0;	c(47)=0;	c(48)=0;	c(49)=0;	c(50)=0;	c(51)=0;	c(52)=0;	c(53)=0;	c(54)=0;
    c(55)=0;	c(56)=0;	c(57)=0;	c(58)=0;	c(59)=0;	c(60)=0;	c(61)=0;	c(62)=0;	c(63)=0;	c(64)=0;	c(65)=0;
    
    PARAMS(2) = psf_pupil;   % size of pupil in mm for which PSF and MTF is to be calculated
    PARAMS(3) = zernike_pupil;	% size of pupil in mm that Zernike coefficients define. NOTE: You can define the aberrations
    % for any pupil size and calculate their effects for any smaller aperture.
    PARAMS(4) = 60*PARAMS(1)*(180/3.1416)*PARAMS(5)*.001/field_size;
    param4orig = PARAMS(4); 	% size of pupil field in mm (use a large field to magnify the PSF)
    
elseif diff_limited == 0;
    [fname,pname] = uigetfile('*.zer','Open Coefficient File');
    fid = fopen([pname fname],'r');
    version = fgetl(fid);% fscanf(fid,'%s',1);
    instrument = fgetl(fid);% fscanf(fid,'%s',1);
    manuf = fgetl(fid);% fscanf(fid,'%s',1);
    oper = fgetl(fid);% fscanf(fid,'%s',1);
    pupoff = fgetl(fid);% fscanf(fid,'%s',1);
    geooff = fgetl(fid);% fscanf(fid,'%s',1);
    datatype = fgetl(fid);% fscanf(fid,'%s',1);
    Rfit = fscanf(fid,'%s %f\n',[1 2]);% fscanf(fid,'%s',1);
    Rfit=Rfit(length(Rfit));
    Rmax = fgetl(fid);% fscanf(fid,'%s',1);
    waverms = fgetl(fid);% fscanf(fid,'%s',1);
    order = fgetl(fid);% fscanf(fid,'%s',1);
    strehl = fgetl(fid);% fscanf(fid,'%s',1);
    refent = fgetl(fid);% fscanf(fid,'%s',1);
    refcor = fgetl(fid);% fscanf(fid,'%s',1);
    resspec = fgetl(fid);% fscanf(fid,'%s',1);
    data = fgetl(fid);% fscanf(fid,'%s',1);
    c = fscanf(fid,'%i %i %g',[3 inf]); %read the first line
    fclose(fid);
    c=c(3,2:66); %ignore the piston term (check this carefully!!!!!)
    
    %     %set these parameter based on what is contained in the *.zer file
    %     PARAMS(2)=2*Rfit; %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
    %     PARAMS(3)=2*Rfit;
    %     PARAMS(4)=15; param4orig = PARAMS(4); %automatically compute the field size
    PARAMS(2) = psf_pupil;   % size of pupil in mm for which PSF and MTF is to be calculated
    PARAMS(3) = zernike_pupil;	% size of pupil in mm that Zernike coefficients define. NOTE: You can define the aberrations
    % for any pupil size and calculate their effects for any smaller aperture.
    PARAMS(4) = 60*PARAMS(1)*(180/3.1416)*PARAMS(5)*.001/field_size;
    param4orig = PARAMS(4); 	% size of pupil field in mm (use a large field to magnify the PSF)
    
else
    %do nothing
end

if (PARAMS(2) < PARAMS(3))
    c = TransformC(c,PARAMS(3),PARAMS(2),0,0,0);
    c(1:5) = 0; %set tilt, defocus and astigmatism to zero
    PARAMS(3) = PARAMS(2); %set PARAMS(3) to correspond to the new coefficients.
end

DefocusInDiopters = defocus;

PARAMS(2)=psf_pupil;

if (PARAMS(2) < PARAMS(3))
    c = TransformC(c,PARAMS(3),PARAMS(2),0,0,0);
    c(1:5) = 0; %set tilt, defocus and astigmatism to zero
    PARAMS(3) = PARAMS(2); %set PARAMS(3) to correspond to the new coefficients.
end

c(4)=(1e6/(4*sqrt(3)))*DefocusInDiopters*((PARAMS(2)/2000)^2); % convert DefocusInDiopters into a Zernike coefficient

% calculate the RMS
rms=sqrt(sum(c(1:65).^2));

% print the result to the screen
%     fprintf('%g\t%g\t',DefocusInDiopters, rms);

%Generate PSF using compPSFOSA.m
PARAMS(4) = param4orig; %reset PARAMS(4) to its original value just in case it was changed in another script
[PSF] = compPSFOSA(c); % call the script compPSFOSA.m to generate PSF

function C2=TransformC(C1,dial,dia2,tx,ty,thetaR)
% “TransformC” returns transformed Zernike coefficient set, C2, from the original set, C1,
% both in standard ANSI order, with the pupil diameter in mm as the first term.
% dia1-pupilk diameter the the coefficient define
% dia2—new pupil diameter [mm]
% tx, ty—Cartesian translation coordinates [mm]
% thetaR—angle of rotation [degrees]
% Scaling and translation is performed first and then rotation.

etaS=dia2/dial; % Scaling factor
etaT=2*sqrt(tx^2+ty^2) / dial; % Translation coordinates
thetaT=atan2(ty, tx);
thetaR=thetaR*pi/180; % Rotation in radians
jnm=length(C1)-1; nmax=ceil((-3+sqrt(9+8*jnm))/2); jmax=nmax*(nmax+3)/2;
S=zeros(jmax+1,1); S(1:length(C1))=C1; C1=S; clear S
P=zeros(jmax+1); % Matrix P transforms from standard to Campbell order 
N=zeros(jmax+1); % Matrix N contains the normalization coefficients
R=zeros(jmax+1); % Matrix R is the coefficients of the radial polynomials
CC1=zeros(jmax+1,1); % CC1 is a complex representation of C1
counter=1;

for m=-nmax:nmax % Meridional indexes
    for n=abs(m):2:nmax % Radial indexes
        jnm=(m+n*(n+2))/2;
        P(counter,jnm+1)=1;
        N(counter,counter)=sqrt(n+1);
        for s=0:(n-abs(m))/2
            R(counter-s,counter)=(-1)^s*factorial(n-s) / (factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s));
        end
    if m<0, CC1(jnm+1)=(C1((-m+n*(n+2))/2+1)+i*C1(jnm+1)) /sqrt(2);
    elseif m==0, CC1(jnm+1)=C1(jnm+1);
    else, CC1(jnm+1)=(C1(jnm+1)-i*C1((-m+n*(n+2))/2+1)) /sqrt(2) ;end
    counter=counter+1;
end, end

ETA=[]; % Coordinate-transfer matrix
for m=-nmax:nmax
    for n=abs(m):2:nmax
        ETA=[ETA P*(transform(n,m,jmax,etaS,etaT,thetaT,thetaR))];
end, end

C=inv(P)*inv(N)*inv(R)*ETA*R*N*P;
CC2=C*CC1;
C2=zeros(jmax+1,1); % C2 is formed from the complex Zernike coefficients, CC2
for m=-nmax:nmax
    for n=abs(m):2:nmax
        jnm=(m+n*(n+2))/2;
        if m<0, C2(jnm+1)=imag(CC2(jnm+1)-CC2((-m+n*(n+2))/2+1)) /sqrt(2);
    elseif m==0, C2(jnm+1)=real(CC2(jnm+1));
    else, C2(jnm+1)=real(CC2(jnm+1)+CC2((-m+n*(n+2))/2+1)) /sqrt(2);
end, end, end
C2=[C2];
%

function Eta=transform(n,m,jmax,etaS,etaT,thetaT,thetaR)
% Returns coefficients for transforming a ro?n*expi*m*theta-term into ’-terms
Eta=zeros(jmax+1,1);
for p=0:((n+m)/2)
    for q=0:((n-m)/2)
        nnew=n-p-q; mnew=m-p+q;
        jnm=(mnew+nnew*(nnew+2))/2;
        Eta(floor(jnm+1))=Eta(floor(jnm+1))+nchoosek((n+m)/2,p)*nchoosek((n-m)/2,q)...
            *etaS^(n-p-q)*etaT^(p+q)*exp(i*((p-q)*(thetaT-thetaR)+m*thetaR));
end, end

function [PSF] = compPSFOSA(c)

global PARAMS;
% Zphase_Mahajan generates the complex pupil function for the Fourier transform

pupilfunc = Zphase_MahajanOSA(c);

pupilfunc=transpose(pupilfunc);

% The amplitude of the point spread function is the Fourier transform
% of the wavefront aberration (when is it expressed in terms of phase)
Hamp=fft2(pupilfunc);

% The intensity of the PSF is the square of the amplitude
% the complex conjugate is a way of multiplying out the complex part of the function

Hint=(Hamp .* conj(Hamp));

% Define the size of the PSF plot in arcmin.
% NOTE: The dimension of a single pixel in the PSF in radians is the wavelength
% divided by the size of the pupil field.
%plotdimension = 60*PARAMS(1)*(180*60/3.1416)*PARAMS(5)*.001/PARAMS(4);

%fprintf('The size of the point spread image is %g seconds of arc\n',plotdimension);

PSF = real(fftshift(Hint)); % this comment reorients the PSF so the origin is at the center of the image
PSF = PSF./(PARAMS(6)^2); % scale the PSF so that peak represents the Strehl ratio
clear Hint;

function phasemap=Zphase_MahajanOSA(c)

% WARNING: 	The Zernike polynomials are normalized and ordered after Mahajan.
%           Applied Optics 33: 8121-8124 (1994)
%				Be careful that your coefficients match the list of polynomials that you use.
%				The terms are ordered so that the sine expressions are always odd.

global PARAMS;

phasemap = zeros(PARAMS(1),PARAMS(1));

sizeoffield=PARAMS(4);

PARAMS(6) = 0;
%*********************************************************************************************
%******************** Parameters for Stiles_Crawford reflectance (From Burns et al)***********
B=0.2;%fraction if the diffuse component;
A=1-B;%fraction of directed component;
peakx=0;%normalized location of peak reflectance in x-direction
peaky=0;%normalized location of peak reflectance in x-direction
p=0;%0.047%set to zero to turn off the amplitude factor (Values from Burns et al in mm^(-2)
%*********************************************************************************************

for ny = 1:PARAMS(1)

	for nx = 1:PARAMS(1)

		xpos = ((nx-1)*(sizeoffield/PARAMS(1))-(sizeoffield/2));
		ypos = ((ny-1)*(sizeoffield/PARAMS(1))-(sizeoffield/2));

      [angle norm_radius]=cart2pol(xpos,ypos);
      norm_radius=norm_radius/(PARAMS(3)/2);
      r=norm_radius;
		%norm_radius = (sqrt(xpos^2+ypos^2))/(PARAMS(3)/2);

		%if (ypos==0 & xpos>0)
		%	angle = 3.1416/2;
		%elseif(ypos==0 & xpos<0)
		%	angle = -3.1416/2;
		%elseif(xpos==0 & ypos==0)
		%	angle = 0;
		%elseif(ypos>0)
		%	angle = atan(xpos/ypos);
		%else
		%	angle= 3.1416 + atan(xpos/ypos);
		%end

	if norm_radius > PARAMS(2)/PARAMS(3)
			waveabermap(nx,ny)=NaN;
	else
			phase = 0;
			phase = ...   
            c(1)*sqrt(4)*((1)*r^1)*sin(1*angle) + ...
            c(2)*sqrt(4)*((1)*r^1)*cos(1*angle) + ...
            c(3)*sqrt(6)*((1)*r^2)*sin(2*angle) + ...
            c(4)*sqrt(3)*((2)*r^2+(-1)*r^0) + ...
            c(5)*sqrt(6)*((1)*r^2)*cos(2*angle) + ...
            c(6)*sqrt(8)*((1)*r^3)*sin(3*angle) + ...
            c(7)*sqrt(8)*((3)*r^3+(-2)*r^1)*sin(1*angle) + ...
            c(8)*sqrt(8)*((3)*r^3+(-2)*r^1)*cos(1*angle) + ...
            c(9)*sqrt(8)*((1)*r^3)*cos(3*angle) + ...
            c(10)*sqrt(10)*((1)*r^4)*sin(4*angle) + ...
            c(11)*sqrt(10)*((4)*r^4+(-3)*r^2)*sin(2*angle) + ...
            c(12)*sqrt(5)*((6)*r^4+(-6)*r^2+(1)*r^0) + ...
            c(13)*sqrt(10)*((4)*r^4+(-3)*r^2)*cos(2*angle) + ...
            c(14)*sqrt(10)*((1)*r^4)*cos(4*angle) + ...
            c(15)*sqrt(12)*((1)*r^5)*sin(5*angle) + ...
            c(16)*sqrt(12)*((5)*r^5+(-4)*r^3)*sin(3*angle) + ...
            c(17)*sqrt(12)*((10)*r^5+(-12)*r^3+(3)*r^1)*sin(1*angle) + ...
            c(18)*sqrt(12)*((10)*r^5+(-12)*r^3+(3)*r^1)*cos(1*angle) + ...
            c(19)*sqrt(12)*((5)*r^5+(-4)*r^3)*cos(3*angle) + ...
            c(20)*sqrt(12)*((1)*r^5)*cos(5*angle) + ...
            c(21)*sqrt(14)*((1)*r^6)*sin(6*angle) + ...
            c(22)*sqrt(14)*((6)*r^6+(-5)*r^4)*sin(4*angle) + ...
            c(23)*sqrt(14)*((15)*r^6+(-20)*r^4+(6)*r^2)*sin(2*angle) + ...
            c(24)*sqrt(7)*((20)*r^6+(-30)*r^4+(12)*r^2+(-1)*r^0) + ...
            c(25)*sqrt(14)*((15)*r^6+(-20)*r^4+(6)*r^2)*cos(2*angle) + ...
            c(26)*sqrt(14)*((6)*r^6+(-5)*r^4)*cos(4*angle) + ...
            c(27)*sqrt(14)*((1)*r^6)*cos(6*angle) + ...
            c(28)*sqrt(16)*((1)*r^7)*sin(7*angle) + ...
            c(29)*sqrt(16)*((7)*r^7+(-6)*r^5)*sin(5*angle) + ...
            c(30)*sqrt(16)*((21)*r^7+(-30)*r^5+(10)*r^3)*sin(3*angle) + ...
            c(31)*sqrt(16)*((35)*r^7+(-60)*r^5+(30)*r^3+(-4)*r^1)*sin(1*angle) + ...
            c(32)*sqrt(16)*((35)*r^7+(-60)*r^5+(30)*r^3+(-4)*r^1)*cos(1*angle) + ...
            c(33)*sqrt(16)*((21)*r^7+(-30)*r^5+(10)*r^3)*cos(3*angle) + ...
            c(34)*sqrt(16)*((7)*r^7+(-6)*r^5)*cos(5*angle) + ...
            c(35)*sqrt(16)*((1)*r^7)*cos(7*angle) + ...
            c(36)*sqrt(18)*((1)*r^8)*sin(8*angle) + ...
            c(37)*sqrt(18)*((8)*r^8+(-7)*r^6)*sin(6*angle) + ...
            c(38)*sqrt(18)*((28)*r^8+(-42)*r^6+(15)*r^4)*sin(4*angle) + ...
            c(39)*sqrt(18)*((56)*r^8+(-105)*r^6+(60)*r^4+(-10)*r^2)*sin(2*angle) + ...
            c(40)*sqrt(9)*((70)*r^8+(-140)*r^6+(90)*r^4+(-20)*r^2+(1)*r^0) + ...
            c(41)*sqrt(18)*((56)*r^8+(-105)*r^6+(60)*r^4+(-10)*r^2)*cos(2*angle) + ...
            c(42)*sqrt(18)*((28)*r^8+(-42)*r^6+(15)*r^4)*cos(4*angle) + ...
            c(43)*sqrt(18)*((8)*r^8+(-7)*r^6)*cos(6*angle) + ...
            c(44)*sqrt(18)*((1)*r^8)*cos(8*angle) + ...
            c(45)*sqrt(20)*((1)*r^9)*sin(9*angle) + ...
            c(46)*sqrt(20)*((9)*r^9+(-8)*r^7)*sin(7*angle) + ...
            c(47)*sqrt(20)*((36)*r^9+(-56)*r^7+(21)*r^5)*sin(5*angle) + ...
            c(48)*sqrt(20)*((84)*r^9+(-168)*r^7+(105)*r^5+(-20)*r^3)*sin(3*angle) + ...
            c(49)*sqrt(20)*((126)*r^9+(-280)*r^7+(210)*r^5+(-60)*r^3+(5)*r^1)*sin(1*angle) + ...
            c(50)*sqrt(20)*((126)*r^9+(-280)*r^7+(210)*r^5+(-60)*r^3+(5)*r^1)*cos(1*angle) + ...
            c(51)*sqrt(20)*((84)*r^9+(-168)*r^7+(105)*r^5+(-20)*r^3)*cos(3*angle) + ...
            c(52)*sqrt(20)*((36)*r^9+(-56)*r^7+(21)*r^5)*cos(5*angle) + ...
            c(53)*sqrt(20)*((9)*r^9+(-8)*r^7)*cos(7*angle) + ...
            c(54)*sqrt(20)*((1)*r^9)*cos(9*angle) + ...
            c(55)*sqrt(22)*((1)*r^10)*sin(10*angle) + ...
            c(56)*sqrt(22)*((10)*r^10+(-9)*r^8)*sin(8*angle) + ...
            c(57)*sqrt(22)*((45)*r^10+(-72)*r^8+(28)*r^6)*sin(6*angle) + ...
            c(58)*sqrt(22)*((120)*r^10+(-252)*r^8+(168)*r^6+(-35)*r^4)*sin(4*angle) + ...
            c(59)*sqrt(22)*((210)*r^10+(-504)*r^8+(420)*r^6+(-140)*r^4+(15)*r^2)*sin(2*angle) + ...
            c(60)*sqrt(11)*((252)*r^10+(-630)*r^8+(560)*r^6+(-210)*r^4+(30)*r^2+(-1)*r^0) + ...
            c(61)*sqrt(22)*((210)*r^10+(-504)*r^8+(420)*r^6+(-140)*r^4+(15)*r^2)*cos(2*angle) + ...
            c(62)*sqrt(22)*((120)*r^10+(-252)*r^8+(168)*r^6+(-35)*r^4)*cos(4*angle) + ...
            c(63)*sqrt(22)*((45)*r^10+(-72)*r^8+(28)*r^6)*cos(6*angle) + ...
            c(64)*sqrt(22)*((10)*r^10+(-9)*r^8)*cos(8*angle) + ...
            c(65)*sqrt(22)*((1)*r^10)*cos(10*angle);

      
         d=sqrt((peakx-xpos)^2+(peaky-ypos)^2);
         SCfactor=B+A*10^(-p*(PARAMS(2)^2)*d^2);
      
      	phasemap(nx,ny) = SCfactor * exp(-i * 2 * 3.1416 * phase/PARAMS(5));
			PARAMS(6) = PARAMS(6) + 1;
		end
	end
end
