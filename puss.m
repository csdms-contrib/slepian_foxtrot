function varargout=puss
% [c,C]=PUSS
%
% OUTPUT:
%
% c    The puss contour
% C    The puss image
% 
% SEE ALSO:
%
% DUCK
%
% Last modified by fjsimons-at-alum.mit.edu, 08/10/2022

% Where I keep the file
fname=fullfile(getenv('EPS'),'puss.mat');

if exist(fname,'file')
  load(fname)
else
  % Get the image
  C=imread(fullfile(getenv('EPS'),'puss.jpg'));
  
  % Black and white
  C=double(round(C/255));
  
  % Contour
  c=contourc(C,[1 1]);
  c=c(:,2:end)';

  % Save
  save(fullfile(getenv('EPS'),'puss.mat'),'c','C')
end
  
% Optional output 
if nargout==0
  imagesc(C)
  colormap(gray(2))
  hold on 
  plot(c(:,1),c(:,2),'y','LineWidth',2)
  hold off
end

varns={c,C};
varargout=varns(1:nargout);



