% THIS PROGRAM DEMONSTRATES HODGKIN HUXLEY MODEL IN CURRENT CLAMP EXPERIMENTS AND SHOWS ACTION POTENTIAL PROPAGATION
% Time is in secs, voltage in mV, conductances in m mho/mm^2, capacitance in uF/mm^2

% Defining the range of impulse currents
irange = 0:0.0001:0.6; 
numPeaksHist = zeros(length(irange)); 

% Defining constants
gkmax = 0.36; 
vk = -77;
gnamax = 1.20;  
vna = 50; 
gl = 0.003; 
vl = -54.387;
cm = 0.01;

% Defining iteration values
dt = 0.01;
niter = 50000;
t = (1:niter) * dt;

% Iterating over the range of currents defined 
for i = 1:length(irange)

    v = -64.9964;
    m = 0.0530;
    h = 0.5960;
    n = 0.3177;
    
    gnahist = zeros(1, niter);
    gkhist = zeros(1, niter);
    vhist = zeros(1, niter);
    mhist = zeros(1, niter);
    hhist = zeros(1, niter);
    nhist = zeros(1, niter);
    
    iapp = irange(i) * ones(1, niter);
    
    % Running the model for each value of current
    for iter = 1:niter
        gna = gnamax * m^3 * h;
        gk = gkmax * n^4;
        gtot = gna + gk + gl;
        
        vinf = ((gna * vna + gk * vk + gl * vl) + iapp(iter)) / gtot;
        tauv = cm / gtot;
        v = vinf + (v - vinf) * exp(-dt / tauv);
        
        alpham = 0.1 * (v + 40) / (1 - exp(-(v + 40) / 10));
        betam = 4 * exp(-0.0556 * (v + 65));
        alphan = 0.01 * (v + 55) / (1 - exp(-(v + 55) / 10));
        betan = 0.125 * exp(-(v + 65) / 80);
        alphah = 0.07 * exp(-0.05 * (v + 65));
        betah = 1 / (1 + exp(-0.1 * (v + 35)));
        
        taum = 1 / (alpham + betam);
        tauh = 1 / (alphah + betah);
        taun = 1 / (alphan + betan);
        
        minf = alpham * taum;
        hinf = alphah * tauh;
        ninf = alphan * taun;
        
        m = minf + (m - minf) * exp(-dt / taum);
        h = hinf + (h - hinf) * exp(-dt / tauh);
        n = ninf + (n - ninf) * exp(-dt / taun);
        
        vhist(iter) = v;
        mhist(iter) = m;
        hhist(iter) = h;
        nhist(iter) = n;
    end
    
    % Counting number of spikes or APs 
    numPeaksHist(i) = spikecounter(vhist);
end

% Finding the indices of positions of spikes 
[ind_i1, ind_i2, ind_i3] = peakcurrentsind(numPeaksHist);
i1 = irange(ind_i1);
i2 = irange(ind_i2);
i3 = irange(ind_i3);
disp("The cutoff currents are:")
disp(i1);
disp(i2);
disp(i3);

% % Plotting characteristics of threshold currents i1, i2, i3

    % figure(1) 
    % plot(t, vhist)
    % title('Voltage vs. Time')
    % xlabel('t')
    % ylabel('V')

    
    % figure(2) 
    % plot(t, mhist, 'y-', t, hhist, 'g.', t, nhist, 'b-')
    % title('Gating Variables vs. Time')
    % legend('m', 'h', 'n')
    % xlabel('t')
    % ylabel('units')

    
    % figure(3)
    % gna = gnamax * (mhist.^3) .* hhist;
    % gk = gkmax * nhist.^4;
    % plot(t, gna, 'r', t, gk, 'b')
    % legend('gna', 'gk')
    % title('Conductance vs. Time')
    % xlabel('t')
    % ylabel('g')

% figure(4)
plot(irange,numPeaksHist*1000/(niter/100));
title('Firing rate vs External current')
xlabel('I(ext)');
ylabel('Firing rate (f)')
hold on;
xline(i1, 'r', 'LineWidth', 2)
hold on;
xline(i2, 'r', 'LineWidth', 2)
hold on;
xline(i3, 'r', 'LineWidth', 2)


% Function defined to count the number of spikes or Action potentials for a
% given current
function numspikes = spikecounter(vhist)
    vthreshold = 10;
    numspikes = 0;
    for i = 2:length(vhist) - 1
        if (vhist(i) > vhist(i - 1)) && (vhist(i) > vhist(i + 1)) && (vhist(i) >= vthreshold)
            numspikes = numspikes + 1;
        end
    end
end

% Function defined to find the index of currents of I1, I2, I3 
function [indi1, indi2, indi3] = peakcurrentsind(numpeaksHist)
    indi1 = 0;
    indi2 = 0;
    indi3 = 0;
    
    for ind = 2:length(numpeaksHist) - 1
        if (numpeaksHist(ind) > 0) && (numpeaksHist(ind - 1) == 0)
            indi1 = ind;
        end
        if (numpeaksHist(ind + 1) - numpeaksHist(ind)) > 4
            indi2 = ind;
        end
        if (numpeaksHist(ind + 1) - numpeaksHist(ind)) < -2
            indi3 = ind;
        end
    end
end


