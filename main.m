clear all
% --------------------Choose landscape--------------------

% LANDSCAPE INDEX:
% Simple Landscape - 1
% Complex Landscape - 2
% Convex Landscape - 3
% Ripple Landscape - 4
% Uneven-Convex Landscape - 5
% 6-Hump Camel Landscape - 6
% Modified DeJong Landscape - 7

landscape = 4;

%% ----------------Plots the chosen landscape----------------
colormap(jet); 
switch landscape
    case 1
        axLimits = [-2 2; -2 2];
        ezmesh(@SimpleLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 40;
    case 2
        axLimits = [-3 7; -3 7];
        ezmesh(@ComplexLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 110;
    case 3
        axLimits = [-10 10; -10 10];
        ezmesh(@ConvexLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 200;
    case 4
        axLimits = [-7 7; -7 7];
        ezmesh(@RippleLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 140;
    case 5
        axLimits = [-10 10; -10 10];
        ezmesh(@UnevenLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 200;
    case 6
        axLimits = [-2 2; -1 1];
        ezmesh(@CamelLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 30;
    case 7
        axLimits = [-2 2; -2 2];
        ezmesh(@ModifiedDejongLandscape,axLimits(1,:),axLimits(2,:))
        plotSize = 40;
end
zlabel("Fitness");
%Makes for more readable code when these are passed as parameters
lowerX = axLimits(1,1);
upperX = axLimits(1,2);
lowerY = axLimits(2,1);
upperY = axLimits(2,2);


%% ----------------------Choose Parameters---------------------
NumSteps=50;
LRate=0.1;
MaxMutate=0.75;
plotFitness = false;
pcolorPlot = false; %Only use w/ init_mode 2
individualStartPt = false; %Still used for GradAsc
initialize_mode = 3;% 1 = Individual (all for HC)
                    % 2 = Systematic Grid
                    % 3 = Random Population

%% ----------------------Main Execution ------------------------------
if initialize_mode == 1
    %--------Select Start Point for individual algorithm Instance----------
    StartPt = [-1.7,2.7];
    %StartPt = randomStartPt(lowerX, upperX);
    disp("Start co-ords: ");
    disp(StartPt);
    %---------------------------------------Select Algorithm--------------------------------------
    %[optimaReached, maxHeight, iterations, points] = GradAscent(StartPt,NumSteps,LRate, landscape, axLimits);
    [optimaReached, maxHeight, iterations, points] = HillClimb(StartPt,NumSteps,MaxMutate, landscape, axLimits);
    drawPath(points, landscape);

elseif initialize_mode == 2
    %----------------------Alternative: Systematic grid of starting positions----------------------
        allPoints = [];
        i = 1;
        gridDensity = 1;
        %1 - Starting point on each integer co-ordinate
        %2 - Starting point on every 0.5 co-oridnate
        for x = gridDensity*lowerX:gridDensity*upperX
            for y = gridDensity*lowerY:gridDensity*upperY
                StartPt = [x/gridDensity y/gridDensity];
                %-----------------------------Select Algorithm-------------------------------------
                %[optimaReached, maxHeight, iterations, points] = GradAscent(StartPt,NumSteps,LRate, landscape, axLimits);
                [optimaReached, maxHeight, iterations, points] = HillClimb(StartPt,NumSteps,MaxMutate, landscape, axLimits);
                if ~ pcolorPlot drawPath(points, landscape), end 
                allPoints = cat(1, allPoints, points);%Used for pcolor plot
                performance(i, 1) = optimaReached;
                performance(i, 2) = iterations;
                i = i + 1;
            end
        end
elseif initialize_mode == 3
    %------------------------------------Random Population-----------------------------------------
    popSize = 50;
    fitnesses = [];
    cumulativeBest = [];
    for x = 1:popSize
        StartPt = randomStartPt(lowerX, upperX);
        [optimaReached, maxHeight, iterations, points] = HillClimb(StartPt,NumSteps,MaxMutate, landscape, axLimits);
        drawPath(points, landscape);
        fitnesses = [fitnesses; maxHeight];
        performance(x, 1) = optimaReached;
        performance(x, 2) = iterations;
    end
    for x = 1:length(fitnesses)
        cumulativeBest(x) = max(fitnesses(:,x));
    end
    figure;
    semilogy(cumulativeBest,'LineWidth',2);
    xlabel('Iteration');
    ylabel('Best Fitness');
    grid on;
end

%Plots best fitness vs iteration
if plotFitness
    figure;
    %plot(BestCost,'LineWidth',2);
    semilogy(maxHeight,'LineWidth',2);
    xlabel('Iteration');
    ylabel('Best Fitness');
    grid on;
end



%% -----------------------Code to generate colourmap plots----------------------------
if pcolorPlot && ~individualStartPt
    %if simpleLandscape plotSize = 40, else plotSize = 110, end
        colourMap = zeros(plotSize);
        %Map tally of points on to coordinate grid for colourmap
        for i = 1:size(allPoints,1)
            %round and scale coordinates appropriately
            x = round((allPoints(i,1)+ abs(lowerX) +.1)*10);
            if x > plotSize x=plotSize; end
            y = round((allPoints(i,2)+ abs(lowerY) +.1)*10);
            if y > plotSize y=plotSize; end
            colourMap(x,y) = colourMap(x,y) + 1;
        end
        %pcolor omits last row and col, add one of each to negate
        newcol = zeros(plotSize,1);
        newrow = zeros(1,plotSize+1);
        colourMap = [colourMap newcol];
        colourMap = [colourMap; newrow];
        c = pcolor(linspace(lowerX,upperX,plotSize+1),linspace(lowerY,upperY,plotSize+1),colourMap');
        c.FaceColor = 'interp';
        colorbar;
end


%% -----------------------Display Appropriate Output------------------------------
if initialize_mode == 1
    disp("Global Optimum Reached: " + optimaReached);
    disp("Maximum Height: " + maxHeight(end));
    disp("Iterations: " + iterations);

else 
    [successProbability, avgIterations] = evaluateSuccess(performance);
    disp("Individual Success Probability: " + successProbability);
    disp("Avg iterations to reach global optima: " +  avgIterations);
end



%% -----------------Practical Function Definitions-----------------------

%A function that draws a line showing the path of chosen algorithm
%Params: coords - 2d array of x,y coords of all points on a path
%        simpleLandscape - Boolean indicating landscape in use
function drawPath(coords, landscape)
    for x = 1:size(coords,1)
        coords(x,3) = calculateHeight(coords(x,1),coords(x,2), landscape);
    end
    hold on
    plot3(coords(:,1), coords(:,2), coords(:,3)+0.1, 'r');
end


%Params: optima - 2d array mapping boolean success to iterations needed
%Returns: - Probability of reaching global optima
%         - Average iterations required to do so
function [successProbability, avgIterations] = evaluateSuccess(optima)
    successes = 0;
    iterations = 0;
    for x = 1:size(optima,1)
        if optima(x,1) == 1
            successes = successes + 1;
            iterations = iterations + optima(x,2);
        end
    end
    successProbability = successes/size(optima,1);
    avgIterations = iterations/successes;
end


%Returns: startPt - Random start point within bounds 
function [startPt] = randomStartPt(lowerBound, upperBound) 
    startPt = lowerBound + rand(1,2)*(upperBound-lowerBound);
end

%Params: x,y - coordinates
%        lansdscape - int representing fitness landscape in use
%Returns: height - the height (fitness) of landscape at x,y
function [height] = calculateHeight(x, y, landscape)
    switch landscape
        case 1
            height = SimpleLandscape(x, y);
        case 2
            height = ComplexLandscape(x, y);
        case 3
            height = ConvexLandscape(x, y);
        case 4
            height = RippleLandscape(x, y);
        case 5
            height = UnevenLandscape(x, y);
        case 6
            height = CamelLandscape(x, y);
        case 7
            height = ModifiedDejongLandscape(x, y);
    end
end


%Function implementing gradient ascent
%Returns: MaxReached - On simple landscape, binary if global optima found
%                    - On complex landscape, max height reached as double
%         Iterations - For simple landscape, iterations to reach G optima
function [optimaReached, maxHeight, iterations, points] = GradAscent(StartPt,NumSteps,LRate, landscape, limits)
    points = [];
	PauseFlag=0;
	hold on;
    optimaReached = 0;
    maxHeight = [0];
    
    
	for i = 1:NumSteps
        
	    %---------Calculates the 'height' at StartPt-----------------
        points = cat(1, points, StartPt);
        height = calculateHeight(StartPt(1),StartPt(2), landscape);
        disp("Height: " + height);
        %Updates max height
        if height > maxHeight(end)
            maxHeight = [maxHeight, height];
        else%To allow plotting of MaxHeight against time (iterations)
            maxHeight = [maxHeight, maxHeight(end)];
        end
        
	    %-----------------Plots point on the landscape-----------------
        plot3(StartPt(1), StartPt(2), height, '*r', 'MarkerSize',10);
        
	
	    %------------Calculates the gradient at StartPt-----------------
        
        switch landscape
            case 1
                grad = SimpleLandscapeGrad(StartPt(1),StartPt(2));
            case 2
                grad = ComplexLandscapeGrad(StartPt(1),StartPt(2));
            case 3
                disp("Derivative missing")
            case 4
                 grad = RippleLandscapeGrad(StartPt(1),StartPt(2));
            case 5
                disp("Derivative missing")
            case 6
                grad = CamelLandscapeGrad(StartPt(1),StartPt(2));
            case 7
                grad = ModifiedDejongLandscapeGrad(StartPt(1),StartPt(2));
        end
        
        
	    %--------Calculates the new point and update StartPt-----------
        StartPt(1) = StartPt(1) + LRate*grad(1);
        StartPt(2) = StartPt(2) + LRate*grad(2);
        
            
	    %---------Ensures StartPt is within the specified bounds--------
        xmin = limits(1,1);
        xmax = limits(1,2);
        ymin = limits(2,1);
        ymax = limits(2,2);
        StartPt = max([StartPt;xmin ymin]);
        StartPt = min([StartPt;xmax ymax]);
            
        %--------------Checks if global optimum reached----------------
        %Overwrites maxReached to binary val if simple landscape
        iterations = i;
        if GlobalOptimaReached(landscape, height)%Checks if global optima reached
            disp("Global Optimum Reached")
            optimaReached = 1;
            break
            %If step limit hit, return boolean 0 for maxima not met
        elseif i == NumSteps %%&& simpleLandscape
            optimaReached = 0;
            break
%         elseif i == NumSteps && ~simpleLandscape
%             maxReached = i;
        end
        
        
	    %-------------------Pauses to view output----------------------
% 	    if(PauseFlag)
% 	        j=input('Press return to continue\nor 0 and return to stop pausing\n');
% 	        if(j==0) PauseFlag=0; end;
%         end
    end
    
    
	hold off
end

% Mutation function
% Params: OldPt - x,y location vector describing current position
%         MaxMutate - the maximum amount by which the current pos can be
%         mutated
% Returns: NewPt- The new, mutated potential location point 
function[NewPt] = Mutate(OldPt,MaxMutate)
	%Selects a random MutDist in the range(-MaxMutate,MaxMutate)
    MutDist = (2*MaxMutate)*rand(1,1) - MaxMutate;
	%---------Mutates random element of OldPt by MutDist-----------
    randElement = randi([1 2],1,1);
    OldPt(randElement) = OldPt(randElement)+ MutDist;
    NewPt = OldPt;
end	

%Returns: MaxReached - On simple landscape, binary if global optima found
%                    - On complex landscape, max height reached as double
%         Iterations - For simple landscape, iterations to reach G optima
function [optimaReached, maxHeight, iterations, points] =  HillClimb(StartPt,NumSteps,MaxMutate, landscape, limits)
    points = [];
	PauseFlag=0;
	hold on;
    optimaReached = 0;
    maxHeight = [0];
    noImprovementCount =0;
    
	for i = 1:NumSteps
        
        %---------Calculates the 'height' at StartPt-----------------
        points = cat(1, points, StartPt);%Comment out for improved pcolor plot
        height = calculateHeight(StartPt(1),StartPt(2), landscape);
        disp("Height: " + height);
        %Updates max height
        if height > maxHeight(end)
            maxHeight = [maxHeight, height];
        else%To allow plotting of MaxHeight against time (iterations)
            maxHeight = [maxHeight, maxHeight(end)];
        end
        
        %---------------Plots point on landscape---------------------
        plot3(StartPt(1), StartPt(2), height, '*r', 'MarkerSize',10);
        
	    %-------------Mutates StartPt in to NewPt-------------------
	    NewPt=Mutate(StartPt,MaxMutate);

        %------Ensures NewPt is within the specified bounds---------
        xmin = limits(1,1);
        xmax = limits(1,2);
        ymin = limits(2,1);
        ymax = limits(2,2);
        NewPt = max([NewPt;xmin ymin]);
        NewPt = min([NewPt;xmax ymax]);

        %------Calculates the height of the new point--------------
        newHeight = calculateHeight(NewPt(1),NewPt(2), landscape);
        disp("New Height: " + newHeight);
                
        %-------Decides whether to update StartPt or not-----------   
        if newHeight >= height
            %points = cat(1, points, StartPt); %un-comment for improved pcolor plot
            StartPt = NewPt;
            disp("Point updated");
        else
            delta = height-newHeight;
            if rand(1,1) > exp(delta/height*4)
                StartPt = NewPt;
                disp("Point updated with less fitness");
            end
            noImprovementCount = noImprovementCount + 1;
        end
        
        %----------------Terminating Conditions--------------------
        %Overwrites maxReached to binary val if simple landscape
        iterations = i;
        if GlobalOptimaReached(landscape, height) %Checks for global optimum
            disp("Global Optimum Reached")
            optimaReached = 1;
            %break %Comment out for Random Pop
        elseif i == NumSteps %%&& simpleLandscape %Checks if iteration limit hit
            optimaReached = 0;
            break
        elseif noImprovementCount == 50000
            optimaReached = 0;
            disp("Stuck in local Minima");
            break
        end
        
	    %-----------------Pauses to view output-----------------
% 	    if(PauseFlag)
% 	        x=input('Press return to continue\nor 0 and return to stop pausing\n');
% 	        if(x==0) PauseFlag=0; end;
% 	    end		
	end
	hold off
end

%Params: landsdcape - int representing fitness landscape in use
%        height - height or fitness of current location
%Return: complete - boolean if sufficient 'global optima' found
function [complete] = GlobalOptimaReached(landscape, height)
    complete = false;
    percentage = 0.95; %Adjust as needed
    switch landscape
            case 1
                globalOptima = percentage*4;
            case 2
                globalOptima = percentage*12.2;
            case 3
                globalOptima = percentage*31.2;
            case 4
                 globalOptima = percentage*8;
            case 5
                globalOptima = percentage*15.6;
            case 6
                globalOptima = percentage*5.8;
            case 7
                globalOptima = percentage*8.1;
    end
    %Checks if current height passes threshold
    if height >= globalOptima
        complete = true;
    end
end



%% -------------------Fitness Landscape Definitions-----------------------

function [z] = SimpleLandscape(x,y)
	z=max(1-abs(2*x),0)+x+y;
end

function [z] = ConvexLandscape(x,y)
	z=((y^2)/16 - (x^2)/25)*5;
end

function [z] = UnevenLandscape(x,y)
	z=cos((x+y))+(x^2)/6-(y^2)/6;
end

function [z] = RippleLandscape(x,y)%Global = 8
	z=(4*cos((x^2+y^2)/4)/(x^2+y^2+1))*2;
end

function [z] = CamelLandscape(x,y)
    term1 = (4-2.1*x^2+(x^4)/3) * x^2;
    term2 = x*y;
    term3 = (-4+4*y^2) * y^2;

    z = term1 + term2 + term3;
end

function [z] = ModifiedDejongLandscape(x,y)
    z = -(x^2 + y^2)+8;
end

function [z] = Stairs(x,y)
	z=(sign(-.65-x)+sign(-.35-x)+sign(-.05-x)+sign(.25-x)+sign(.55-x))/7;
end

function [f]=ComplexLandscape(x,y) %Global Optima: 12.2039
	f=4*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*(x/5 - x^3 - y^5)*exp(-x^2-y^2) -(1/3)*exp(-(x+1)^2 - y^2)-1*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9)*exp(-(x-3)^2-(y-3)^2);
end

function z=Rosenbrock(x)
    n=numel(x);
    z=sum((1-x(1:n-1)).^2)+100*sum((x(2:n)-x(1:n-1).^2).^2);
end

function [y] = rastr(x,y)
sum = 0;

sum = sum + (x^2 - 10*cos(2*pi*x));


y = 10*2 + sum;

end

%% ---------------------Landscape Gradient Functions-----------------------

function [g] = ModifiedDejongLandscapeGrad(x, y)
    g(1) = -2*x;
    g(2) = -2*y;
end

function z=Rastrigin(x)
    
    n=numel(x);

    A=10;
    
    z=n*A+sum(x.^2-A*cos(2*pi*x));

end

function [g] = CamelLandscapeGrad(x,y)
    g(1) = 2*(x^5-(4.2)*x^3+4*x+0.5*y);
    g(2) = x+16*y^3-8*y;
end

function [g] = RippleLandscapeGrad(x,y)
    g(1) = (-4*x*((x^2+y^2+1)*sin((1/4)*(x^2+y^2))+4*cos((1/4)*(x^2+y^2))))/(x^2+y^2+1)^2;
    g(2) = (-4*y*((x^2+y^2+1)*sin((1/4)*(x^2+y^2))+4*cos((1/4)*(x^2+y^2))))/(x^2+y^2+1)^2;
end

function [g] = SimpleLandscapeGrad(x,y) %Global Optima: 4
	if(1-abs(2*x) > 0)
	    if(x<0) g(1) = 3;
	    elseif(x==0) g(1)=0;
	    else g(1) = -1;
	    end
	else g(1)=1;
	end
	g(2)=1;
end

% Definition of gradient of Complex landscape
function [g]=ComplexLandscapeGrad(x,y)
	g(1)=-8*exp(-(x^2)-(y+1)^2)*((1-x)+x*(1-x)^2)-15*exp(-x^2-y^2)*((0.2-3*x^2) -2*x*(x/5 - x^3 - y^5)) +(2/3)*(x+1)*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*(14*(x-3)^6-2*(x-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
	g(2)=-8*(y+1)*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*exp(-x^2-y^2)*(-5*y^4 -2*y*(x/5 - x^3 - y^5)) +(2/3)*y*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*((-1.5*(y-4)^4+9*(y-3)^8)-2*(y-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
end