function [num, den, z, p] = butter(n, Wn, varargin)
%BUTTER Butterworth digital and analog filter design.
%   [B,A] = BUTTER(N,Wn) designs an Nth order lowpass digital
%   Butterworth filter and returns the filter coefficients in length 
%   N+1 vectors B (numerator) and A (denominator). The coefficients 
%   are listed in descending powers of z. The cut-off frequency 
%   Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to 
%   half the sample rate.
%
%   If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an 
%   order 2N bandpass filter with passband  W1 < W < W2.
%   [B,A] = BUTTER(N,Wn,'high') designs a highpass filter.
%   [B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
%   
%   When used with three left-hand arguments, as in
%   [Z,P,K] = BUTTER(...), the zeros and poles are returned in
%   length N column vectors Z and P, and the gain in scalar K. 
%
%   When used with four left-hand arguments, as in
%   [A,B,C,D] = BUTTER(...), state-space matrices are returned.
%
%   BUTTER(N,Wn,'s'), BUTTER(N,Wn,'high','s') and BUTTER(N,Wn,'stop','s')
%   design analog Butterworth filters.  In this case, Wn can be bigger
%   than 1.0.
%
%   See also BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ, FILTER.

%   Author(s): J.N. Little, 1-14-87
%   	   J.N. Little, 1-14-88, revised
%   	   L. Shure, 4-29-88, revised
%   	   T. Krauss, 3-24-93, revised

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

[btype,analog,errStr] = iirchk(Wn,varargin{:});
error(errStr)

if n>500
	error('Filter order too large.')
end

% step 1: get analog, pre-warped frequencies
if ~analog,
	fs = 2;
	u = 2*fs*tan(pi*Wn/fs);
else
	u = Wn;
end

Bw=[];
% step 2: convert to low-pass prototype estimate
if btype == 1	% lowpass
	Wn = u;
elseif btype == 2	% bandpass
	Bw = u(2) - u(1);
	Wn = sqrt(u(1)*u(2));	% center frequency
elseif btype == 3	% highpass
	Wn = u;
elseif btype == 4	% bandstop
	Bw = u(2) - u(1);
	Wn = sqrt(u(1)*u(2));	% center frequency
end

% step 3: Get N-th order Butterworth analog lowpass prototype
[z,p,k] = buttap(n);

% Transform to state-space
[a,b,c,d] = zp2ss(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
if btype == 1		% Lowpass
	[a,b,c,d] = lp2lp(a,b,c,d,Wn);

elseif btype == 2	% Bandpass
	[a,b,c,d] = lp2bp(a,b,c,d,Wn,Bw);

elseif btype == 3	% Highpass
	[a,b,c,d] = lp2hp(a,b,c,d,Wn);

elseif btype == 4	% Bandstop
	[a,b,c,d] = lp2bs(a,b,c,d,Wn,Bw);
end

% step 5: Use Bilinear transformation to find discrete equivalent:
if ~analog,
	[a,b,c,d] = bilinear(a,b,c,d,fs);
end

if nargout == 4
	num = a;
	den = b;
	z = c;
	p = d;
else	% nargout <= 3
% Transform to zero-pole-gain and polynomial forms:
	if nargout == 3
		[z,p,k] = ss2zp(a,b,c,d,1);
		z = buttzeros(btype,n,Wn,Bw,analog);
		num = z;
		den = p;
		z = k;
	else % nargout <= 2
		den = poly(a);
		num = buttnum(btype,n,Wn,Bw,analog,den);
		% num = poly(a-b*c)+(d-1)*den;

	end
end

%---------------------------------
function b = buttnum(btype,n,Wn,Bw,analog,den)
% This internal function returns more exact numerator vectors
% for the num/den case.
% Wn input is two element band edge vector
if analog
    switch btype
    case 1  % lowpass
        b = [zeros(1,n) n^(-n)];
        b = real( b*polyval(den,-j*0)/polyval(b,-j*0) );
    case 2  % bandpass
        b = [zeros(1,n) Bw^n zeros(1,n)];
        b = real( b*polyval(den,-j*Wn)/polyval(b,-j*Wn) );
    case 3  % highpass
        b = [1 zeros(1,n)];
        b = real( b*den(1)/b(1) );
    case 4  % bandstop
        r = j*Wn*((-1).^(0:2*n-1)');
        b = poly(r);
        b = real( b*polyval(den,-j*0)/polyval(b,-j*0) );
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
    case 1  % lowpass
        r = -ones(n,1);
        w = 0;
    case 2  % bandpass
        r = [ones(n,1); -ones(n,1)];
        w = Wn;
    case 3  % highpass
        r = ones(n,1);
        w = pi;
    case 4  % bandstop
        r = exp(j*Wn*( (-1).^(0:2*n-1)' ));
        w = 0;
    end
    b = poly(r);
    % now normalize so |H(w)| == 1:
    kern = exp(-j*w*(0:length(b)-1));
    b = real(b*(kern*den(:))/(kern*b(:)));
end

function z = buttzeros(btype,n,Wn,Bw,analog)
% This internal function returns more exact zeros.
% Wn input is two element band edge vector
if analog
    % for lowpass and bandpass, don't include zeros at +Inf or -Inf
    switch btype
    case 1  % lowpass
        z = zeros(0,1);
    case 2  % bandpass
        z = zeros(n,1);
    case 3  % highpass
        z = zeros(n,1);
    case 4  % bandstop
        z = j*Wn*((-1).^(0:2*n-1)');
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
    case 1  % lowpass
        z = -ones(n,1);
    case 2  % bandpass
        z = [ones(n,1); -ones(n,1)];
    case 3  % highpass
        z = ones(n,1);
    case 4  % bandstop
        z = exp(j*Wn*( (-1).^(0:2*n-1)' ));
    end
end

% FUNCTIONS

function [at,bt,ct,dt] = lp2lp(a,b,c,d,wo)
%LP2LP Lowpass to lowpass analog filter transformation.
%   [NUMT,DENT] = LP2LP(NUM,DEN,Wo) transforms the lowpass filter
%   prototype NUM(s)/DEN(s) with unity cutoff frequency of 1 rad/sec 
%   to a lowpass filter with cutoff frequency Wo (rad/sec).
%   [AT,BT,CT,DT] = LP2LP(A,B,C,D,Wo) does the same when the
%   filter is described in state-space form.
%
%   See also BILINEAR, IMPINVAR, LP2BP, LP2BS and LP2HP

%   Author(s): J.N. Little and G.F. Franklin, 8-4-87

if nargin == 3		% Transfer function case
        % handle column vector inputs: convert to rows
        if size(a,2) == 1
            a = a(:).';
        end
        if size(b,2) == 1
            b = b(:).';
        end
	% Transform to state-space
	wo = c;
	[a,b,c,d] = tf2ss(a,b);
end

error(abcdchk(a,b,c,d));
[ma,nb] = size(b);
[mc,ma] = size(c);

% Transform lowpass to lowpass
at = wo*a;
bt = wo*b;
ct = c;
dt = d;

if nargin == 3		% Transfer function case
    % Transform back to transfer function
    [z,k] = tzero(at,bt,ct,dt);
    num = k * poly(z);
    den = poly(at);
    at = num;
    bt = den;
end


function [at,bt,ct,dt] = lp2bp(a,b,c,d,wo,bw)
%LP2BP	Lowpass to bandpass analog filter transformation.
%	[NUMT,DENT] = LP2BP(NUM,DEN,Wo,Bw) transforms the lowpass filter
%	prototype NUM(s)/DEN(s) with unity cutoff frequency to a
%	bandpass filter with center frequency Wo and bandwidth Bw.
%	[AT,BT,CT,DT] = LP2BP(A,B,C,D,Wo,Bw) does the same when the
%	filter is described in state-space form.


if nargin == 4		% Transfer function case
	% Transform to state-space
	wo = c;
	bw = d;
	[a,b,c,d] = tf2ss(a,b);
end

error(abcdchk(a,b,c,d));
[ma,nb] = size(b);
[mc,ma] = size(c);

% Transform lowpass to bandpass
q = wo/bw;
at = wo*[a/q eye(ma); -eye(ma) zeros(ma)];
bt = wo*[b/q; zeros(ma,nb)];
ct = [c zeros(mc,ma)];
dt = d;

if nargin == 4		% Transfer function case
	% Transform back to transfer function
	b = poly(at);
	at = poly(at-bt*ct)+(dt-1)*b;
	bt = b;
end


function [at,bt,ct,dt] = lp2hp(a,b,c,d,wo)
%LP2HP	Lowpass to highpass analog filter transformation.
%	[NUMT,DENT] = LP2HP(NUM,DEN,Wo) transforms the lowpass filter
%	prototype NUM(s)/DEN(s) with unity cutoff frequency to a
%	highpass filter with cutoff frequency Wo.
%	[AT,BT,CT,DT] = LP2HP(A,B,C,D,Wo) does the same when the
%	filter is described in state-space form.

if nargin == 3		% Transfer function case
	% Transform to state-space
	wo = c;
	[a,b,c,d] = tf2ss(a,b);
end

error(abcdchk(a,b,c,d));
[ma,nb] = size(b);
[mc,ma] = size(c);

% Transform lowpass to highpass
at =  wo*inv(a);
bt = -wo*(a\b);
ct = c/a;
dt = d - c/a*b;

if nargin == 3		% Transfer function case
	% Transform back to transfer function
	b = poly(at);
	at = poly(at-bt*ct)+(dt-1)*b;
	bt = b;
end


function [at,bt,ct,dt] = lp2bs(a,b,c,d,wo,bw)
%LP2BS	Lowpass to bandstop analog filter transformation.
%	[NUMT,DENT] = LP2BS(NUM,DEN,Wo,Bw) transforms the lowpass filter
%	prototype NUM(s)/DEN(s) with unity cutoff frequency to a
%	bandstop filter with center frequency Wo and bandwidth Bw.
%	[AT,BT,CT,DT] = LP2BS(A,B,C,D,Wo,Bw) does the same when the
%	filter is described in state-space form.

if nargin == 4		% Transfer function case
	% Transform to state-space
	wo = c;
	bw = d;
	[a,b,c,d] = tf2ss(a,b);
end

error(abcdchk(a,b,c,d));
[ma,nb] = size(b);
[mc,ma] = size(c);

% Transform lowpass to bandstop
q = wo/bw;
at =  [wo/q*inv(a) wo*eye(ma); -wo*eye(ma) zeros(ma)];
bt = -[wo/q*(a\b); zeros(ma,nb)];
ct = [c/a zeros(mc,ma)];
dt = d - c/a*b;

if nargin == 4		% Transfer function case
	% Transform back to transfer function
	b = poly(at);
	at = poly(at-bt*ct)+(dt-1)*b;
	bt = b;
end


function [zd, pd, kd, dd] = bilinear(z, p, k, fs, fp, fp1)
%BILINEAR Bilinear transformation with optional frequency prewarping.
%   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
%   function specified by Z, P, and K to a z-transform discrete
%   equivalent obtained from the bilinear transformation:
%
%      H(z) = H(s) |
%                  | s = 2*Fs*(z-1)/(z+1)
%
%   where column vectors Z and P specify the zeros and poles, scalar
%   K specifies the gain, and Fs is the sample frequency in Hz.
%   [NUMd,DENd] = BILINEAR(NUM,DEN,Fs), where NUM and DEN are 
%   row vectors containing numerator and denominator transfer
%   function coefficients, NUM(s)/DEN(s), in descending powers of
%   s, transforms to z-transform coefficients NUMd(z)/DENd(z).
%   [Ad,Bd,Cd,Dd] = BILINEAR(A,B,C,D,Fs) is a state-space version.
%   Each of the above three forms of BILINEAR accepts an optional
%   additional input argument that specifies prewarping. For example,
%   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs,Fp) applies prewarping before
%   the bilinear transformation so that the frequency responses
%   before and after mapping match exactly at frequency point Fp
%   (match point Fp is specified in Hz).
%
%   See also IMPINVAR.

%   Author(s): J.N. Little, 4-28-87 
%   	   J.N. Little, 5-5-87, revised

%   Gene Franklin, Stanford Univ., motivated the state-space
%   approach to the bilinear transformation.

[mn,nn] = size(z);
[md,nd] = size(p);

if (nd == 1 & nn < 2) & nargout ~= 4	% In zero-pole-gain form
	if mn > md
		error('Numerator cannot be higher order than denominator.')
	end
	if nargin == 5		% Prewarp
		fp = 2*pi*fp;
		fs = fp/tan(fp/fs/2);
	else
		fs = 2*fs;
	end
	z = z(finite(z));	 % Strip infinities from zeros
	pd = (1+p/fs)./(1-p/fs); % Do bilinear transformation
	zd = (1+z/fs)./(1-z/fs);
% real(kd) or just kd?
	kd = (k*prod(fs-z)./prod(fs-p));
	zd = [zd;-ones(length(pd)-length(zd),1)];  % Add extra zeros at -1

elseif (md == 1 & mn == 1) | nargout == 4 %
	if nargout == 4		% State-space case
		a = z; b = p; c = k; d = fs; fs = fp;
		error(abcdchk(a,b,c,d));
		if nargin == 6			% Prewarp
			fp = fp1;		% Decode arguments
			fp = 2*pi*fp;
			fs = fp/tan(fp/fs/2)/2;
		end
	else			% Transfer function case
		if nn > nd
			error('Numerator cannot be higher order than denominator.')
		end
		num = z; den = p;		% Decode arguments
		if nargin == 4			% Prewarp
			fp = fs; fs = k;	% Decode arguments
			fp = 2*pi*fp;
			fs = fp/tan(fp/fs/2)/2;
		else
			fs = k;			% Decode arguments
		end
		% Put num(s)/den(s) in state-space canonical form.  
		[a,b,c,d] = tf2ss(num,den);
	end
	% Now do state-space version of bilinear transformation:
	t = 1/fs;
	r = sqrt(t);
	t1 = eye(size(a)) + a*t/2;
	t2 = eye(size(a)) - a*t/2;
	ad = t2\t1;
	bd = t/r*(t2\b);
	cd = r*c/t2;
	dd = c/t2*b*t/2 + d;
	if nargout == 4
		zd = ad; pd = bd; kd = cd;
	else
		% Convert back to transfer function form:
		p = poly(ad);
		zd = poly(ad-bd*cd)+(dd-1)*p;
		pd = p;
	end
else
	error('First two arguments must have the same orientation.')
end

