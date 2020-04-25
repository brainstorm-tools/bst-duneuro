%%
function [condtensor, s, fail] = sb_calcTensorCond_tuch_wm(cfg,mask,V1,V2,V3,L1,L2,L3)
L1.img=1000*L1.anatomy;%turn it to um^2/msec from mm^2/sec
L2.img=1000*L2.anatomy;
L3.img=1000*L3.anatomy;

index_wm = 1;
index_csf = 2;

check1 = isequal(V1.dim,V2.dim,V3.dim);
check2 = isequal(V1.dim,size(V1.anatomy),size(V2.anatomy),size(V3.anatomy));
check3 = isequal(L1.dim,L2.dim,L3.dim);
check4 = isequal(L1.dim,size(L1.anatomy),size(L2.anatomy),size(L3.anatomy));
check5 = isequal(V1.dim(1:3),L1.dim,mask.dim);
check6 = ~(isempty(cfg)||isempty(cfg.conductivity));
% check7 = ~(length(cfg.conductivity)<5);
check7 = 1;
if (check1 && check2 && check3 && check4 && check5 && check6 && check7)
    fail       = 0;
    failnan    = 0;
    condtensor = cell(mask.dim);
    N = zeros(1,3);
    vol = zeros(1,3);
    for i = 1 : V1.dim(1)
        for j = 1 : V1.dim(2)
            for k = 1 : V1.dim(3)
                if ((mask.anatomy(i,j,k) == index_wm))&(L1.anatomy(i,j,k)>10^-4)&(L2.anatomy(i,j,k)>10^-4)&(L3.anatomy(i,j,k)>10^-5) % the L values are significants >> then zeros
                    S = zeros(3);
                    S(:,1) = V1.anatomy(i,j,k,:); % get the three vectors .... main direction
                    S(:,2) = V2.anatomy(i,j,k,:);
                    S(:,3) = V3.anatomy(i,j,k,:);
                    if(norm(S'*S-diag([1,1,1]),2)>10e-7) %% ++++++++++ In this case the computation fails the vecors are not orthogonal
                        S(:,1) = S(:,1) / norm(S(:,1));
                        S(:,2) = S(:,2) - (S(:,1)'*S(:,2))*S(:,1);
                        S(:,2) = S(:,2) / norm(S(:,2));
                        S(:,3) = S(:,3) - (S(:,1)'*S(:,3))*S(:,1) - (S(:,2)'*S(:,3))*S(:,2);
                        S(:,3) = S(:,3) / norm(S(:,3));
                        failnan = failnan + 1;
                    end
                    if(sum(sum(isnan(S),1),2)>0) % the orientation is not well defined ... the tensor will be remain isotrop and no translation
                        condtensor{i,j,k} = cfg.conductivity(mask.anatomy(i,j,k))*diag([1,1,1]);
                    else % +++++ This case seems to be the isotrope case .... the tensor is isotropic and then translate to the right coordinate system  ++++ %
                        D = diag([L1.anatomy(i,j,k),L2.anatomy(i,j,k),L3.anatomy(i,j,k)]); 
                        T = S * D * S';
                        %T = hier skalieren
%                         T
                        condtensor{i,j,k} = T;
                        %                         if (mask.anatomy(i,j,k) == 5)
                        %                             N(1) = N(1) + 1;
                        %                             vol(1) = vol(1) + L1.anatomy(i,j,k)*L2.anatomy(i,j,k)*L3.anatomy(i,j,k);
                        if (mask.anatomy(i,j,k) == index_wm) % ++++ get the WM volume and compute 
                            N(1) = N(1) + 1;
                            vol(1) = vol(1) + L1.anatomy(i,j,k)*L2.anatomy(i,j,k)*L3.anatomy(i,j,k);
                            %                    elseif (mask.anatomy(i,j,k) == 3)
                            %                        N(3) = N(3) + 1;
                            %                        vol(3) = vol(3) + L1.anatomy(i,j,k)*L2.anatomy(i,j,k)*L3.anatomy(i,j,k);
                        end
                    end
                elseif(((mask.anatomy(i,j,k) == index_wm)))
                    condtensor{i,j,k} = cfg.conductivity(mask.anatomy(i,j,k))*diag([1,1,1]);
                    if(( mask.anatomy(i,j,k) == index_wm )&(~((L1.anatomy(i,j,k)>10^-4) && ...
                            (L2.anatomy(i,j,k)>10^-4)&(L3.anatomy(i,j,k)>10^-5))))
                        fail = fail + 1;
                    end
                else % +++++ This case seems to be the air where the conductivity will be set to zero  ++++ %
                    condtensor{i,j,k} = zeros(3);
                end
            end
        end
    end
    
    d(1) = vol(1) / N(1); % volume moyen des tissu anisotrop
    d(1) = d(1)^(1/3); % cubic root ... equivalent to a distance now... mean distance
    %d(2) = vol(2) / N(2);
    %d(2) = d(2)^(1/3);
    %d(3) = vol(3) / N(3);
    %d(3) = d(3)^(1/3);
    s = d(1)*cfg.conductivity(index_wm);%+d(2)*cfg.conductivity(6); ===> conductivity absolute ... pas S/mm mais just Siemens 
    s = s / (d(1)^2);% + d(2)^2)); %  =====> densite de sonductivity moyenne par surface  S/m.m ==> later it will be multiplied by factor of distance tenosr
    fprintf('s*d = %.6f\n',s*d(1))
    fprintf('s = %.6f\n',s)
    %    failper = fail / (N(1)+N(2));
    mx = -10000;
    mn = 10000;
    for i = 1 : V1.dim(1)
        for j = 1 : V1.dim(2)
            for k = 1 : V1.dim(3)
                if ((mask.anatomy(i,j,k) == 1)&(L1.anatomy(i,j,k)>10^-4)&(L2.anatomy(i,j,k)>10^-4)&(L3.anatomy(i,j,k)>10^-5))
                    condtensor{i,j,k} = s * condtensor{i,j,k};
                    % +++++++++ Adapt the diffusion tensor to the conductivity tensor by linesr multiplication by the mean factor value of variation mappd to the conductivity value of the white matter  
                    %keep outliers wma bigger from the highest cond around cond of wm
                    if max(condtensor{i,j,k}(:)) > cfg.conductivity(index_csf)
                        condtensor{i,j,k}(1:1+size(condtensor{i,j,k},1):end) = cfg.conductivity(index_wm);
                    end
                    
                    if mx < max(max(condtensor{i,j,k}))
                        mx = max(max(condtensor{i,j,k}));
                        %fprintf('mx %d, %d, %d\n',i,j,k)
                    end
                    if mn > min(min(condtensor{i,j,k}))
                        mn = min(min(condtensor{i,j,k}));
                        %fprintf('mn %d, %d, %d\n',i,j,k)
                    end
                end
            end
        end
    end
    disp(mn)
    disp(mx)
end
end