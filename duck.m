function varargout=duck(p,t)
% [d,D]=DUCK(p,t)
%
% INPUT:
%
% p    A color threshold [defaulted]
% t    A cutoff row [defaulted]
%
% OUTPUT:
%
% d    The duck's contour
% D    The duck's image
% 
% SEE ALSO:
%
% PUSS
% 
% Last modified by fjsimons-at-alum.mit.edu, 08/09/2022

try 
  load(fullfile(getenv('EPS'),'duck.mat'));
catch
  defval('p',30)
  defval('t',417)

  % Get the image
  D=imread(fullfile(getenv('EPS'),'duck.jpg'));

  % Black and white
  D(D>p)=255;
  D(D<=p)=0;

  % Reflection point
  D(t:end,:,:)=deal(255);

  % Sum it all up
  D=nanmean(double(D),3);

  % Contour
  d=contourc(D,[255 255]);
  d=d(:,2:end)';

  % Save
  save(fullfile(getenv('EPS'),'duck.mat'),'d','D')
end

% Optional output 
if nargout==0
  image(D)
  colormap(gray(2))
  hold on 
  plot(d(:,1),d(:,2),'y','LineWidth',2)
  hold off
end

varns={d,D};
varargout=varns(1:nargout);
