% Show Jupiter's orbit around the sun
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES THAT CAN BE CHANGED
% slices: Number of angles used, this includes 0 and 2pi
% num_of_equal_areas: splits ellipse into equal area sections. << slices.
% Keep num_of_equal_areas much less than slics for best results
% Higher num_of_equal sizes shows better results, still << slices though.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slices = 500;
num_of_equal_areas = 25;                % number of equal areas we want.



% Properties of the ellipse of the planet (Jupiter)
a = 778.570e6;      % semi-major axis (km)
b = 776.874e6;      % semi-minor axis (km)
e=0.9;%e = 0.0489;         % eccentricity
P = 11.86;          % period of jupiter in years
% Per = 740.522e6;    % smallest distnace from sun to jupiter (perihilion)
% Aph = 816.618e6;    % largest distance from sun to jupiter (Aphelion)

% Period to radians
P = P*(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we can set up the orbit using what we know about jupiter and the
% rules of an ellipse. Second we'll try to use kepler's 2nd law to relate
% the position of jupiter to time passed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slices = 500;
theta = linspace(0,2*pi,slices);         % angle set
r = a*(1-e^2)./(1+e*cos(theta));         % distance set
figure(1)
% % % % % % polarplot(theta,r/1e6)                   % Trajectory of pluto around sun.
[x,y] = pol2cart(theta,r/1e6);
% plot(x,y);
% hold on

E = 2*atan(sqrt((1-e)/(1+e)).*tan(theta./2));   % Eccentric anomaly
A = abs(0.5*a*b.*(E-e*sin(E)));                 % Areas from 0 to theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixing A vector so it shows total area traversed from 0 to 2pi for both
% even and odd slices case.

% Odd 'slices' case
iseven = ~bitget(slices,1);     % 1 if even, 0 if odd
if(iseven == 0)
    m = floor(slices/2)+1;      % marks the elementposition for pi.
    A2 = A;                     % for comparison of methods
    for n = m+1:slices
        A(n) = (A(m) - A(n)) + A(m);  
        A2(n) = (pi*a*b) - A2(n);     % correct result probably
    end
end

% Even 'slices' case
if(iseven == 1)
    m = slices/2;      % marks the elementposition for pi.
    mid_slice = pi/2*a*b - A(m);
    half_area = A(m) + mid_slice;
    A2 = A;                     % for comparison of methods
    for n = m+1:slices
        A(n) = (half_area - A(n)) + half_area;   % both work now (A and A2)
        A2(n) = (pi*a*b) - A2(n);                % correct result probably
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second we decide how many equal areas we would like to split the orbit
% into. Using more physics based algorithms would make this unnecessary but
% since we are just using ellipses we'll need to find another way.

% The proposed method is to designate some number of equal areas we want to
% split the orbit into, preferibly the number of equal areas we want should
% be much smaller than the number of slices we have. The smaller the ratio
% the better the results should be. 

% This is because what we will do is check each area element (A) and
% compare to the area needed for equal number of areas. We then stop once
% we start moving away from require area. We use a counter to know how many
% slices we traversed to get near that equal area number. Then this counter
% is used to change how fast we update the plot showing the orbit. So the
% number of equal areas is the number of 'speeds' jupiter will move at and
% this should look like it is accelarating and decelarating. This works
% because of Kepler's 2nd Law, or 3rd Law, where equal areas are traversed
% in equal amount of time. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% num_of_equal_areas = 47;                % number of equal areas we want.
dA_E = pi*a*b / num_of_equal_areas;     % the area of each equal section.

counter = zeros(1,num_of_equal_areas);
diff = zeros(1,2);                      % 1 hold upper diff, 2 holds lower diff
c = 1;
L = 1;

% While loop that counts how many slices are traversed to cover equal area
% partions, this number is held in dA_E.
n = 1;
while n <= slices-1
    upper = A(n+1) - A(L);
    lower = A(n) - A(L);
    if(upper > dA_E && lower <= dA_E)
        diff(1) = abs(dA_E - upper);
        diff(2) = abs(dA_E - lower);
        if(diff(1) > diff(2))                % if upper is farther than lower
            counter(c) = n-L; 
            L = n;
            n = n - 1;                       % this way we start again at the same n.
        else
            counter(c) = n+1-L;
            L = n+1;
        end
        c = c + 1;
    end
    n = n+1;
end
% This is in case the final value isn't reached, we then just set a speed
% because it definitlye shouldn't be 0.
if counter(num_of_equal_areas) == 0
    counter(num_of_equal_areas) = n-L;
end
% This deals with small numbers at the end by making it at least as fast as
% the one before, ideally it would be speeding up. 
if (counter(num_of_equal_areas) < counter(num_of_equal_areas - 1))
    counter(num_of_equal_areas) = counter(num_of_equal_areas - 1);
end

% Error message will play if slice to equal areas ratio causes bad results.
if(length(counter) ~= num_of_equal_areas)
    errmsg = ["counter() should be of length num_of_equal_areas. If counter()" ...,
        "length is higher, then increase slices or decrease num_of_equal_areas."];
    warning(errmsg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third we use what we've gathered to make plots and record a video.
% Using horchler's QTWriter program allows us to output each frame to a
% .mov video file but more importantly allows us to determent the frame
% rate of each frame so we end up with a video with variable frame rate.
% 
% This is important because The method used to get the orbit was not
% physics based but geometric. Instead we find how many slices we must 
% traverse until we cover another dA_E, then these slices crossed tell us
% the relative speed to the other dA_E's crossed so we know how often to 
% update the plot. With QTWriter we are able to just make the framerate be
% equal to the number of slices for that dA_E section. We used Kepler's 
% Law that states equal areas are crossed at equal times. 
% 
% horchler's QTWriter program
% https://github.com/horchler/QTWriter.git
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hf = figure(1);
title("Slices: " + slices + "  Equal Part: " + num_of_equal_areas);
ylabel("10^6 km");
xlabel("10^6 km");
hold on
orbits = 0; 

% Plot in a polar plot. Also uncomment the one in the while loop.
% polarplot(theta,r/1e6,'b');
% polarplot(0,0,'y-o');


% Plot on cartesian coordinates
plot(x,y);
plot(0,0,'y-o');

movObj = QTWriter('orbit.mov');
while orbits < 2 

    % What c and cc do needs to be made clear. Say counter(cc) is 5. Then using
    % c, we can make counter(cc) be applied 5 times. if counter(cc) were 8,
    % then it would be applied 8 times. By this I mean that the framerate of
    % the next 8 frames will be 8. It turns counter = [5,8] the useful form
    % counter = [5,5,5,5,5,8,8,8,8,8,8,8,8].
    c = 1;                
    cc = 1;
    for n = 1:slices-1

%         N = polarplot(theta(n),r(n)/(1e6), 'k-o');
        N = plot(x(n),y(n), 'k-o');
        drawnow;
        if(counter(cc) == c && cc ~= num_of_equal_areas)
            cc = cc+1;
            c = 0;
        end
        
        c = c + 1;
        movObj.FrameRate = counter(cc);
        writeMovie(movObj,getframe(hf));
        delete(N);                            % remove old planet location
    end
    orbits = orbits + 1;
end
close(movObj);

