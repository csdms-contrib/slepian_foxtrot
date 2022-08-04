function varargout=svdslep3(XY,KXY,J,ngro,tol,xver)
% [E,V,c11cmnR,c11cmnK,SE,KXY]=SVDSLEP3(XY,KXY,J,ngro,tol,xver)
%
% Two-dimensional Slepian functions with arbitrary concentration/limiting
% regions in the Cartesian spatial and (half-)spectral domains.
%
% INPUT:
%
% XY       [X(:) Y(:)] coordinates of a SPATIAL-domain curve, in PIXELS
% KXY      [X(:) Y(:)] coordinates of a SPECTRAL-domain HALF-curve, 
%          i.e. in the POSITIVE (including zero) spectral halfplane. 
%          The coordinates here are RELATIVE to the ultimate size of that
%          half-plane, with the CORNER points assuming UNIT wavenumber 
%          values, so they will be ratios, fractions of the Nyquist
%          plane, after the computational growth of the SPACE plane
% J        Number of eigentapers requested [default: 10] 
% ngro     The computational "growth factor" [default: 3]
% tol      abs(log10(tolerance)) for EIGS [default: 12]
% xver     Performs excessive verification, making plots
%
% OUTPUT:
%
% E        The eigenfunctions of the concentration problem
% V        The eigenvalues of the concentration problem
% c11cmnR  The spatial coordinates of the top left corner after growth
% c11cmnK  The spectral coordinates of the top left corner after growth
% SE       The periodogram of the eigenfunctions
% XY       The input spatial-domain curve
% KXY      The symmetrized spectral-space domain curve
%
% EXAMPLE:
%
% svdslep3('demo1',ngro) % with ngro the growth factor
% 
% SEE ALSO:
%
% LOCALIZATION2D
%
% Last modified by fjsimons-at-alum.mit.edu, 08/02/2022

% Default values
defval('J',10);
defval('ngro',4);
defval('tol',12);
defval('xver',1);

% Default SPATIAL curve is a CIRCLE in PIXEL space, of radius cR and cN points
defval('cR',30)
defval('cN',41)
defval('XY',...
       cR*[cos(linspace(0,2*pi,cN)) ; sin(linspace(0,2*pi,cN))]')

% And some (half-)square in the SPECTRAL (half-)space
% Here you use the notion of the Shannon ratio as in SVDSLEP2, which is
% relative to the grown area which has unit (kx,ky) in the CORNERS
defval('R',0.1)
% Remember this R is strictly for convenience in the next line
defval('KXY',...
       R*[-1 1  1 -1 -1; 1 1  0  0  1]')

if ~isstr(XY)
  % Check if we've already computed these
  % Make a hash with the input variables so you don't recompute
  fname=hash([XY(:)' KXY(:)' J ngro tol],'SHA-256');
  % You're going to need an environmental variable and an appropriate directory
  fnams=fullfile(getenv('IFILES'),'HASHES',sprintf('%s_%s.mat',upper(mfilename),fname));

  % Compute and save or load if presaved
  if ~exist(fnams,'file') | 1==1
    t=tic;

    % Check the curves and return the range on the inside 
    % For the SPATIAL part, before the growth domain, in pixels
    [XY,xlimt,ylimt,Nx,Ny]=ccheck(XY,0,[1 1],xver*0);

    % Now embed this in a larger-size matrix to get rid of edge effects and
    % control the spectral discretization. The spatial domain is the unique
    % reference, in pixels, for everything that comes below. Mind the order!
    newsize=ngro*[Ny Nx];
    
    % Expand the SPATIAL domain and recompute the coordinates
    [QinR,c11cmnR,QX,QY]=cexpand(XY,xlimt,ylimt,[Ny Nx],newsize,0,xver*0);

    % For the SPECTRAL part, mirrored, in the discretization appropriate to
    % the growth domain, still relative to the Nyquist plane
    [KXY,kxlimt,kylimt,Nkx,Nky]=ccheck(KXY,1,1./newsize,xver*0);

    % Expand the SPECTRAL domain to return the full Nyquist plane
    [QinK,c11cmnK,QKX,QKY]=cexpand(KXY,kxlimt,kylimt,[Ny Nx],newsize,1,xver*0);

    % Enforce Hermiticity by taking out the first row and column of match
    % since the dci component (see KNUM2) will be in the lower right
    % quadrant since the data set will be even no matter what
    [i,j]=ind2sub(newsize,QinK);
    QinK=QinK(i>min(i)&j>min(j));
    
    % The result must be Hermitian!
    dom=zeros(newsize); dom(QinK)=1;
    difer(isreal(ifft2(ifftshift(dom)))-1,[],[],[])

    % Now you are ready for a check
    if xver==1
      clf
      % Plot the SPATIAL region of interest exactly as input
      ah(1)=subplot(221);
      twoplot(XY,'b'); hold on
      % Compute the centroid - my function is from SAGA, but R2022a has its own
      [X0,Y0]=centroid(XY(:,1),XY(:,2));
      plot(X0,Y0,'b+')
      axis equal; grid on
      % Just open up the axes a bit
      xlim(minmax(XY(:))*1.5)
      ylim(minmax(XY(:))*1.5)

      % Plot the SPECTRAL region of interest after the mirroring operation
      ah(2)=subplot(222);
      twoplot([KXY ; KXY(1,:)],'r'); hold on
      % This centroid remains zero
      [KX0,KY0]=deal(0);
      plot(KX0,KY0,'r+')
      axis equal; grid on
      % Open up the axes to the full Nyquist plane, corners at [+/-1 +/-1]
      axis([-1 1 -1 1])

      ah(3)=subplot(223);
      % Plot the SPATIAL region of interest after the growth domain
      plot(QX(QinR),QY(QinR),'k.')
      % The original curve in PIXEL coordinates needs to plot right on here
      hold on
      twoplot(XY,'k','LineWidth',2)
      hold off
      axis equal; grid on

      ah(4)=subplot(224);
      % Plot the SPECTRAL region of interest after the growth domain
      plot(QKX(QinK),QKY(QinK),'k.')
      % The curve in FRACTIONAL coordinates needs to plot right on here
      hold on
      % Fiddle with the curve by one unit
      twoplot([KXY ; KXY(1,:)],'b','LineWidth',2)
      hold off
      grid on
      axis([-1 1 -1 1])

      ntix=4;
      % Non-pathological Pixel coordinates should round gracefully
      set(ah(3),'xtick',sort(unique(round(...
       		 [c11cmnR(1):round(range(c11cmnR([1 3]))/ntix):c11cmnR(3) c11cmnR(3)]))))
      set(ah(3),'ytick',sort(unique(round(...
       		 [c11cmnR(4):round(range(c11cmnR([2 4]))/ntix):c11cmnR(2) c11cmnR(2)]))))
      % Spectral coordinates are already fractions and need to include% zero
      set(ah(4),'xtick',sort(unique(round(...
       		 [0 c11cmnK(1):round(range(c11cmnK([1 3])*100)/ntix)/100:c11cmnK(3) c11cmnK(3)]*100)/100)))
      set(ah(4),'ytick',sort(unique(round(...
       		 [0 c11cmnK(4):round(range(c11cmnK([2 4]))*100/ntix)/100:c11cmnK(2) c11cmnK(2)]*100)/100)))
      set(ah(3:4),'GridLineStyle',':')

      disp(sprintf('\nHit ENTER to proceed\n'))
      pause
    end

    % The SPECTRAL domain this needs to be turned into a Fourier operator
    QinK=indeks(fftshift(v2s(1:prod(newsize))),QinK);

    % Now make the operators that we are trying to diagonalize
    P=@(x) proj(x,QinR);
    % We're finding VECTORS that are going to be 2-D images!
    Q= @(x) fft2vec(x);
    Qi=@(y) ifft2vec(y);
    L=@(y) proj(y,QinK);
    H=@(x) P(Qi(L(Q(P(x)))));

    % And then find the eigenvectors and eigenvalues
    OPTS.isreal=false;
    OPTS.disp=0;
    defval('tolerance',10^-tol);
    OPTS.tol=tolerance;
    OPTS.maxit=500;

    % Remember to specify the output size
    [E,V]=eigs(H,prod(newsize),J,'LR',OPTS);
    
    [V,i]=sort(diag(V),'descend');
    E=E(:,i); V=V(1:J); E=E(:,1:J);

    % Define some kind of tolerance level
    tol=sqrt(eps); 

    % Make them real as we know they should be
    if any(imag(V)>tol)
      error('Complex eigenvalues');
    else
      V=real(V);
      % Check imaginary part of the "good" eigenfunctions
      disp(sprintf('mean(abs(imag(E))) = %8.3e out of %8.3e\n',...
		   mean(mean(abs(imag(E(:,V>tol))))),...
		   mean(mean(abs(E(:,V>tol))))))
      % Note that they were normalized in the complex plane
      E=real(E); E=E./repmat(diag(sqrt(E'*E))',size(E,1),1);
    end

    if nargout>4
      % Get the power spectrum
      SE=zeros(prod(newsize),size(E,2));
      for i=1:size(E,2),
	SE(:,i)=indeks((abs(fft2(v2s(E(:,i)))).^2),':');
      end
    else
      SE=NaN;
    end

    disp(sprintf('%s took %f seconds',upper(mfilename),toc(t)))
    save(fnams,'E','V','c11cmnR','c11cmnK','SE','XY','KXY')
  else
    disp(sprintf('%s loading %s',upper(mfilename),fnams))
    load(fnams)
  end
  % Output
  varns={E,V,c11cmnR,c11cmnK,SE,XY,KXY};
  varargout=varns(1:nargout);
elseif strcmp(XY,'demo1')
  % Fake the second input now as the growth factor
  defval('KXY',[]); ngro=KXY; clear KXY

  % A circle in SPACE...
  cR=30;
  cN=41;
  XY=cR*[cos(linspace(0,2*pi,cN)) ; sin(linspace(0,2*pi,cN))]';

  % A random blob, fix the radius to be something sizable in pixels
  %[x,y]=blob(1,3); XY=[x y]*20; 

  % And a BOX in SPECTRAL space, no need to close it as it will get
  % mirrored anyway about the lower symmetry axis...
  R=0.1;
  KXY=R*[-1 -1 1 1 ; 0 1 1 0]';

  % A random blob in relative coordinates that are somewhat appropriate
  %[kx,ky]=blob(1,3); KXY=[kx ky]/3; KXY=KXY(KXY(:,2)>=0,:);

  % How many eigenfunctions?
  J=60;
  % Compute the eigenfunctions
  [E,V,c11cmnR,c11cmnK,SE,XY,KXY]=svdslep3(XY,KXY,J,ngro);

  % Plot the first offs+3 basis functions
  offs=0;

  % Make the figures
  figure(1)
  clf
  [ah,ha]=krijetem(subnum(2,3));
  for ind=1:3
    % SPACE-domain functions in PIXEL units
    axes(ah(ind))
    imagefnan(c11cmnR(1:2),c11cmnR(3:4),v2s(E(:,ind+offs)))
    hold on
    plot(XY(:,1),XY(:,2),'k'); hold off
    title(sprintf('%s = %i | %s = %7.5f','\alpha',ind+offs,'\lambda',V(ind+offs)))
    xlabel('horizontal pixels')
    ylabel('vertical pixels')

    % SPECTRAL-domain functions, periodogram
    axes(ha(2*ind))
    psdens=fftshift(decibel(v2s(SE(:,ind+offs))));
    psdens(psdens<-80)=NaN;
    imagefnan(c11cmnK(1:2),c11cmnK(3:4),psdens);
    hold on
    % Remember the original curve was relative to the Nyquist plane
    plot(KXY(:,1),KXY(:,2),'b','LineWidth',2); hold off
    xlabel('scaled horizontal wavenumbers')
    ylabel('scaled vertical wavenumbers')
  end
  fig2print(gcf,'landscape')

  % Also try this one here
  figure(2)
  clf
  EE=sum(repmat(V(:)',length(E),1).*E.^2,2);
  SEE=sum(repmat(V(:)',length(E),1).*SE.^2,2);

  % SPACE-domain functions in PIXEL units
  subplot(121)
  imagefnan(c11cmnR(1:2),c11cmnR(3:4),v2s(EE)); axis image 
  hold on; plot(XY(:,1),XY(:,2),'k'); hold off
  title('Eigenvalue weighted SPATIAL sum')
  xlabel('horizontal pixels')
  ylabel('vertical pixels')

  % SPECTRAL-domain functions, periodogram
  subplot(122)
  psdens=fftshift(decibel(v2s(SEE)));
  psdens(psdens<-80)=NaN;
  imagefnan(c11cmnK(1:2),c11cmnK(3:4),psdens); axis image
  hold on
  % Remember the original curve was relative to the Nyquist plane
  plot(KXY(:,1),KXY(:,2),'k'); hold off
  title('Eigenvalue weighted SPECTRAL sum')
  xlabel('scaled horizontal wavenumbers')
  ylabel('scaled vertical wavenumbers')

  % Also try this one here
  figure(3)
  clf
  plot(V,'o-')
  title(sprintf('sum of the eigenvalues %8.3f',sum(V)))
  longticks(gca,2)
  ylim([-0.1 1.1])
  grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=ccheck(XY,isk,qdxdy,xver)
% Given PIXEL coordinates of a closed curve XY, makes a symmetric centered
% enclosing grid. For the SPATIAL domain (isk=0), that's it (qdxdy=1). For
% the SPECTRAL domain (isk=1), it assumes the curve is in the half-plane,
% and mirrors it, but you now require qdxdy. Stripped from LOCALIZATION2D.
%
% OUTPUT:
%
% XY            The input curve (isk=0), or its mirrored version (isk=1)
% xlimt,ylimt   The limits of the grid
% Nx,Ny         The dimensions of the grid

% In the spectral domain the input curve will be mirrored
defval('isk',0)
% This has been uber-verified as it is
defval('xver',0);

% Make sure the XY of the curve has two columns
if ~isempty(find(size(XY,2)==2))
  if size(XY,1)<size(XY,2)
    XY=XY';
  end
else
  error('Coordinates of the curve not as expected')
end

% Find limits in x and y so as to contain the curves
if isk==0
  % We're in x-space
  xlimt=minmax(XY(:,1));
  ylimt=minmax(XY(:,2));
  % The spacings are in pixels
  qdxdy=[1 1];
elseif isk==1
  % We're in k-space, mirror the half curve
  xlimt=[-1 1]*max(abs(XY(:,1)));
  % Note that input Y must be >= 0
  ylimt=[-1 1]*max(XY(:,2));
  % Actually mirror the curve; do not stick in NaNs or they might end up
  % not matching any input zeros; may need to watch POLY2CW
  XY=[XY ; -XY];
end

% Calculate the full range of the input
xrange=xlimt(2)-xlimt(1);
yrange=ylimt(2)-ylimt(1);

% Grid spacings
qdx=qdxdy(1);
qdy=qdxdy(2);

% You'd use these if you wanted a rectangular grid
Nx=round(xrange/qdx);
Ny=round(yrange/qdy);
% but we make the dimensions (not the domain!) square for now
[Nx,Ny]=deal(max(Nx,Ny));

% Only need to do this if you want to inspect it
if xver==1
  % Don't be a fanatic about half pixels as in LOCALIZATION2D but 
  % simply strive for symmetry
  qx=linspace(xlimt(1),xlimt(2),Nx);
  qy=linspace(ylimt(1),ylimt(2),Ny);
  
  % Remember that these may not be square depending on the choices above
  [QX,QY]=meshgrid(qx,qy);

  % The "curve" is not the boundary but rather the last set of "elements" on the "grid".
  % The midpoint indices of this subset that fall inside of the region...
  Qin=find(inpolygon(QX,QY,XY(:,1),XY(:,2)));
  
  % Make a plot
  figure(1)
  clf
  plot(QX,QY,'k.')
  hold on; axis image
  plot(QX(Qin),QY(Qin),'bo')
  hold off
  xlim(minmax(QX(:))*1.1)
  ylim(minmax(QY(:))*1.1)
  shrink(gca,1.25,1.25); longticks(gca,2)
  t=title(sprintf('Verifying CCHECK %i x %i isk %i',...
		  size(QX),isk));
  movev(t,max(abs(QY(:)))/20)
  disp(sprintf('\nHit ENTER to proceed or CTRL-C to abort\n'))
  % Plot the original curve in actual coordinates
  hold on
  plot(XY(:,1),XY(:,2),'y','LineWidth',2)
  hold off
  pause
end

% Optional output
varns={XY,xlimt,ylimt,Nx,Ny};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qin,c11cmn,QX,QY]=cexpand(XY,xlimt,ylimt,oldsize,newsize,isk,xver)
% Expands a rectangular area enclosing a curve to a new size and computes
% the indices of the interior points and the axis limits

defval('isk',0)
% This also should work like a charm all the time
defval('xver',1);

% How many rows and columns to add on either size?
addon=round([newsize-oldsize]/2);

if isk==0
  % Actually add to the space domain in pixel spacing
  addx=range(ylimt)/oldsize(2)*addon(2);
  addy=range(xlimt)/oldsize(1)*addon(1);
  c11=[xlimt(1) ylimt(2)]+[-addx  addy];
  cmn=[xlimt(2) ylimt(1)]+[ addx -addy];
else
  % The Nyquist plane in wavenumber space in relative coordinates
  c11=[-1  1];
  cmn=[ 1 -1];
end
c11cmn=[c11 cmn];

% Now compute the coordinates in the embedding
qx=linspace(c11(1),cmn(1),newsize(2));
qy=linspace(cmn(2),c11(2),newsize(1));
[QX,QY]=meshgrid(qx,qy);
Qin=find(inpolygon(QX,QY,XY(:,1),XY(:,2)));

if xver==1
  figure(2)
  clf
  plot(QX,QY,'k.')
  hold on; axis image
  plot(QX(Qin),QY(Qin),'bo')
  hold off
  xlim(minmax(QX(:))*1.1)
  ylim(minmax(QY(:))*1.1)
  shrink(gca,1.25,1.25); longticks(gca,2)
  t=title(sprintf('Verifying CEXPAND %i x %i from %i x %i isk %i',...
		  newsize,oldsize,isk));
  movev(t,max(abs(QY(:)))/20)
  disp(sprintf('\nHit ENTER to proceed or CTRL-C to abort\n'))
  % Plot the original curve in actual coordinates
  hold on
  plot(XY(:,1),XY(:,2),'y','LineWidth',2)
  hold off
  pause
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pv=proj(v,p)
% Projects the vector v on the indices p
Pv=zeros(size(v));
Pv(p)=v(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fv=fft2vec(v)
% Returns the two-dimensional FFT of a vector
Fv=fft2(v2s(v));
Fv=Fv(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iFv=ifft2vec(Fv)
% Returns the two-dimensional IFFT of a vector
iFv=ifft2(v2s(Fv));
iFv=iFv(:);

