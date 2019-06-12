clear all
close all

load star_snapshots
load star_data

%% find the target stars as the center of each grid points at t=0;
% inc=180-data(:,3);
r_stars = [x(:,181) y(:,181) z(:,181)];
theta_stars = atan2(r_stars(:,2),r_stars(:,1));
r_stars_norm = vecnorm(r_stars')';
% star_target_grids = zeros(30,32);
% star_database=[];
% 
% for radius = 2:1:31 % Move up by 1000 stars to be settled
% 
%         r_mid= radius+0.5;
%         num_star=round(5000* 2*r_mid / (32^2-2^2));
%         
%         
%     for theta = -pi: 2*pi/num_star :(pi-2*pi/num_star)
%         r_stars = [x(:,181) y(:,181) z(:,181)];
%         
%         cond1=r_stars_norm>=radius;
%         cond2=r_stars_norm<(radius+1);
%         cond3=theta_stars>=theta;
%         cond4=theta_stars<(theta+2*pi/num_star);
%         cond5= inc<10;
%         
%         id = find(cond1 & cond2 & cond3 & cond4 &cond5);
%         % find star closest to grid center
%         
%         theta_mid=theta+pi/num_star;
%         r_grid_center =[r_mid*cos(theta_mid) r_mid*sin(theta_mid) 0];
%         
%         i_all=1:1e5+1;
%         i_bad= setdiff(i_all,id);
%         r_stars(i_bad',:)=repmat([inf inf inf],length(i_bad),1);
%         
%         idx= knnsearch(r_stars,r_grid_center,'K',1);
%         
%         star_database=[star_database idx-1];
%         
%         hold on
%         plot3(x(idx,181),y(idx,181),z(idx,181),'. k','MarkerSize',20)
%         hold on
%         plot3(r_grid_center(1),r_grid_center(2),r_grid_center(3),'*')
%         
%         
%     end
% end
%%

global R_min R_max
R_min=2;
R_max=32;

theta= data(:,end);

load star_target_database_inc_10
j_t=[];

for j=2:32
    
    for k=0:31
        
        idr=find(r_stars_norm>j & r_stars_norm<j+1);
        idth=find(theta>-180+(k)*180/16 & theta<-180+(k+1)*180/16);
        id=intersect(idr,idth);
        id1=intersect(star_database,id-1);
        
        if isempty(id)
            id=0;
        end
        %     id=find(r_stars_norm<k);
%         id1=setdiff(star_database,id-1);
        error_J_term = J_N_r_theta(id1',data);
        j_t=[j_t; k error_J_term];
    end
end

% figure(3)
% plot(j_t(:,1),j_t(:,2))
% xlabel('theta below which stars are removed')
% ylabel('J')


figure(4)
plot(j_t(:,2))
xlabel('theta below which stars are removed')
ylabel('J')
