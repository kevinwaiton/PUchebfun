%1D function on [-1 1]

F1 = PUchebfun(@(x) atan(100*(x-0.8)));


plot(F1); hold on;
show(F1,-1.5,0.65);
axis square; axis tight;
title('atan(100(x-0.8)) with subdomains plotted');

disp('Press any key to continue')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
pause;

close();

%2D function on [0 1;0 1]

F2 = PUchebfun(@(x,y) log(1+10^5*((x-0.5).^2+(y-0.5).^2)),[0 1;0 1]);

subplot(1,2,1);
plot(F2);
axis square; axis tight;
title('log(1+10^5((x-0.5)^2+(y-0.5)^2))');
subplot(1,2,2);
axis square; axis tight;
title('subdomain plots');
show(F2);

disp('Press any key to continue')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
pause;

close();

%3D function on [-1 1;-1 1;-1 1]

F3 = PUchebfun(@(x,y,z) atan(5*(x+y+z)));


