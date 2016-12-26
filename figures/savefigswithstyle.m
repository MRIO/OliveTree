% savefigswithstyle.m



function savefigswithstyle(gcf,'style')

     fnam='your_fig.png'; % your file name

% the engine
% ...get style sheet info

     snam='1col'; % note: NO extension...
     s=hgexport('readstyle',snam);

     s.Format = 'png'

% ...apply style sheet info
     hgexport(gcf,fnam,s);