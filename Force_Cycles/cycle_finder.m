% M+M's 3 thru 10-cycle finder

function [out3,out4,out5,out6,out7,out8,out9,out10] = cycle_finder(varargin)

% Input:
% A = adjacency matrix
% coords = matrix of coordinates for all particles with the center set at 
% (0,0) (optional; only used in practice for 2D simulations)

A = varargin{1};
if nargin > 1
    coords = varargin{2};
    xlength = max(coords(:,1)) - min(coords(:,1));
    ylength = max(coords(:,2)) - min(coords(:,2));
else
    coords = [];
    xlength = 0;
    ylength = 0;
end

count3 = 0;
count4 = 0;
count5 = 0;
count6 = 0;
count7 = 0;
count8 = 0;
count9 = 0;
count10 = 0;
cyc3 = zeros(1e6,3);
cyc4 = zeros(1e6,4);
cyc5 = zeros(1e6,5);
cyc6 = zeros(1e6,6); 
cyc7 = zeros(1e6,7);
cyc8 = zeros(1e6,8);
cyc9 = zeros(1e6,9);
cyc10 = zeros(1e6,10); % Initialize lists of cycles

for i = 1:length(A) % Loop over all reference particles
    if sum(A(i,:)) < 2
        continue; % Skip particles with 0 or 1 network pair
    end
    disp(i);
    ni = find(A(i,:)~=0); % Find all the links to the reference particle
    for j = ni(ni > i) % Loop over all particles linked to r.p.
        if sum(A(j,:)) == 1
            continue; % Check to see if this is the end of a chain
        end
        nj = find(A(j,:)~=0); % Find all the links with this second particle
        for k = nj(nj > i) % Repeat this process as needed; shorten with a while loop?
            if sum(A(k,:)) == 1 % Check to see if this is the end of a chain
                continue; % If so, move on to the next contact
            end
            nk = find(A(k,:)~=0);
            if sum(i == nk) > 0 % Check that the cycle closes on itself
               temp = [i j k];
               if true_cycle(temp,A,coords,xlength,ylength)
                   count3 = count3+1;
                   cyc3(count3,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                   continue; % Move on to the next contact
               end 
            elseif nargout > 1    
                for l = nk(nk > i)
                    if sum(A(l,:)) == 1
                        continue; % Check to see if this is the end of a chain
                    end
                    nl = find(A(l,:)~=0);
                    if sum(i == nl) > 0 % Check that the cycle closes on itself
                        temp = [i j k l];
                        if true_cycle(temp,A,coords,xlength,ylength)
                            count4 = count4+1;
                            cyc4(count4,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                            continue; % Move on to the next contact
                        end
                    elseif nargout > 2    
                        for m = nl(nl > i)
                            if sum(A(m,:)) == 1
                                continue; % Check to see if this is the end of a chain
                            end
                            nm = find(A(m,:)~=0);
                            if sum(i == nm) > 0 % Check that the cycle closes on itself
                                temp = [i j k l m];
                                if true_cycle(temp,A,coords,xlength,ylength)
                                    count5 = count5+1;
                                    cyc5(count5,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                                    continue; % Move on to the next contact
                                end
                            elseif nargout > 3
                                for n = nm(nm > i)
                                    if sum(A(n,:)) == 1
                                        continue; % Check to see if this is the end of a chain
                                    end
                                    nn = find(A(n,:)~=0);
                                    if sum(i == nn) > 0 % Check that the cycle closes on itself
                                        temp = [i j k l m n];
                                        if true_cycle(temp,A,coords,xlength,ylength)
                                            count6 = count6+1;
                                            cyc6(count6,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                                            continue; % Move on to the next contact
                                        end
                                    elseif nargout > 4
                                        for p = nn(nn > i)
                                            if sum(A(p,:)) == 1
                                                continue; % Check to see if this is the end of a chain
                                            end
                                            np = find(A(p,:)~=0);
                                            if sum(i == np) > 0 % Check that the cycle closes on itself
                                                temp = [i j k l m n p];
                                                if true_cycle(temp,A,coords,xlength,ylength) 
                                                    count7 = count7+1;
                                                    cyc7(count7,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                                                    continue; % Move on to the next contact
                                                end
                                            elseif nargout > 5
                                                for q = np(np > i)
                                                    if sum(A(q,:)) == 1
                                                        continue; % Check to see if this is the end of a chain
                                                    end
                                                    nq = find(A(q,:)~=0);
                                                    if sum(i == nq) > 0 % Check that the cycle closes on itself
                                                        temp = [i j k l m n p q];
                                                        if true_cycle(temp,A,coords,xlength,ylength)  
                                                            count8 = count8+1;
                                                            cyc8(count8,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                                                            continue; % Move on to the next contact
                                                        end
                                                    elseif nargout > 6
                                                        for r = nq(nq > i)
                                                            if sum(A(r,:)) == 1
                                                                continue; % Check to see if this is the end of a chain
                                                            end
                                                            nr = find(A(r,:)~=0);
                                                            if sum(i == nr) > 0
                                                                temp = [i j k l m n p q r];
                                                                if true_cycle(temp,A,coords,xlength,ylength)
                                                                    count9 = count9+1;
                                                                    cyc9(count9,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                                                                    continue; % Move on to the next contact
                                                                end
                                                            elseif nargout > 7
                                                                for s = nr(nr > i)
                                                                    if sum(A(s,:)) == 1
                                                                        continue; % Check to see if this is the end of a chain
                                                                    end
                                                                    ns = find(A(s,:)~=0);
                                                                    if sum(i == ns) > 0
                                                                        temp = [i j k l m n p q r s];
                                                                        if true_cycle(temp,A,coords,xlength,ylength)
                                                                            count10 = count10+1;
                                                                            cyc10(count10,:) = sort(temp); % If all conditions are satisfied, add to the list of cycles
                                                                            continue; % Move on to the next contact
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

out3 = unique(cyc3(1:count3,:),'rows');
out4 = unique(cyc4(1:count4,:),'rows');
out5 = unique(cyc5(1:count5,:),'rows');
out6 = unique(cyc6(1:count6,:),'rows'); 
out7 = unique(cyc7(1:count7,:),'rows');
out8 = unique(cyc8(1:count8,:),'rows');
out9 = unique(cyc9(1:count9,:),'rows');
out10 = unique(cyc10(1:count10,:),'rows'); % Eliminate identical cycles
