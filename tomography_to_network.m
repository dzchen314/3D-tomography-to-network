%------------------------------------------------------------------------
%This program imports 3d sliced images and processes them into labeled
%particles and individualized contacts
%Created by David Chen 5/17/17
%email: ender314@gmail.com
%------------------------------------------------------------------------

%Parameters:
%boxheight in meters
boxheight = 0.105;

%ahi is highest image #, alo is lowest image #
ahi = 720;
alo = 1;
numimages = ahi-alo+1; %number of images

filename = sprintf('/filepath/scan_1/slice_50.jpg');%test image to establish image size and voxel dimensions
I = imread(filename);
%resample images down to improve performance and rescale for correct
%isotropic voxel dimensions (zcal = zdim calibration factor for this set)
zcal = numimages/(size(I,2)/24.3*18*numimages/720);%720 is total number of images taken in 180mm
I = imresize(I,[round(zcal*size(I,1)) round(zcal*size(I,2))]);

%looping to process all the scans, s = scan #, a = image #
for s=1:5
    
II=zeros(size(I,1),size(I,2),numimages); %preallocating an array to store the 3d images

    for a=alo:ahi

    %import images from files
    filename = sprintf('/filepath/scan_%d/slice_%d.jpg',s,a);
    I = imread(filename);
    I = imresize(I,[round(zcal*size(I,1)) round(zcal*size(I,2))]);

    %preprocess the images (FFT filtering)
    clear yy;
    clear xx;
    I=medfilt2(I);

    F = fft2(I);
    Fm = fftshift(F);
    Fm = abs(Fm);

    Fp = cos(angle(F)) + 1j*sin(angle(F));

    Fsize = size(F);

    theta = 25*2*pi/360; %25 degree tolerance for lines
    offsetfromcenter = round(Fsize(1,1)/25);
    y1limit = round(Fsize(1,1)/2)-offsetfromcenter;
    y2limit = round(Fsize(1,1)/2)+offsetfromcenter;
    for yy = 1:y1limit
        for xx = round(Fsize(1,2)/2-(y1limit-yy)*atan(theta)):round(Fsize(1,2)/2+(y1limit-yy)*atan(theta))
            if Fm(yy,xx)>10000
                Fm(yy,xx)=0*Fm(yy,xx)/3;
            end
        end
    end
    clear yy;
    clear xx;
    for yy = y2limit:Fsize(1,1)
        for xx = round(Fsize(1,2)/2-(yy-y2limit)*atan(theta)):round(Fsize(1,2)/2+(yy-y2limit)*atan(theta))
            if Fm(yy,xx)>10000
                Fm(yy,xx)=0*Fm(yy,xx)/3;
            end
        end
    end

    Fm = ifftshift(Fm);
    F = Fm.*Fp;
    filtered = mat2gray(abs(ifft2(F)));
    %nlfilt = NLMF(filtered); - optional non-local means filtering. requires a non-local means filter package
    I=filtered;

    %  Options.kernelratio=1;
    %  Options.windowratio=1;
    %  Options.filterstrength=0.1;

    I = imadjust(I);
    I = mat2gray(I);

    I = imbinarize(I,0.4);
    I = medfilt2(I);
    I = ~bwareaopen(~I,50);
    I = bwareaopen(I,50);
    %imtool(I)
    
    %finished image preprocessing

    II(:,:,a) = I;

        if any(a==20:20:500)
            fprintf('loading images, percent complete: %.2f\n',(a)/(500)*100);
        end

    end

%perform distance transform on 3d image stack
fprintf('Calculating distance transform on stack %d (~7 s)\n',s);
tic
IIi = imcomplement(II);
D = bwdist(IIi,'euclidean');
D = -D;
toc

%find regional minima in inverted distance transform and superimpose for
%watershed (marker-based watershed)
fprintf('Calculating and affixing regional minima (~3 min)\n');
tic
Min = imextendedmin(D,5); %5 here needs adjusting depending on how many total pixels you have
LD = imimposemin(D,Min);
toc

%make bright regions infinity to improve watershed calculation
LD(~II) = inf;

%perform watershed
fprintf('Calculating watershed on stack %d (~10 min)\n',s);
tic
W = watershed(LD); %check quality of this before proceeding!!!
toc

%extract labeled particles
fprintf('Extracting and calculating particle properties on stack %d (~<1 min)\n',s);
tic
particles = W;
particles(II==0)=0;

%extract contacts
contacts = II&W==0;

%extract particle properties
cstats = regionprops(contacts,'Area','Centroid');
carea = [cstats.Area];
ccenter = [cstats.Centroid];

%convert particles from solid spheres to hollow surfaces (reduces data)
psurfaces = bwperim(particles);
psurfaces = bwlabeln(psurfaces);

psstats = regionprops(psurfaces,'Area','Centroid','FilledArea');
psarea = [psstats.Area];
parea = [psstats.FilledArea];
pscenter = [psstats.Centroid];

psurf = psurfaces;
toc

%initialize arrays to store particle information
numparts = max(max(max(psurf())));
partcents = zeros(numparts,3);
meshsize = 100;
partlocs = NaN((meshsize+2)*(meshsize*2+2),3,numparts);
totalpartvol = 0;
partscounted = 0;

%generate surface mesh of particles, remove outliers, and perform smoothing
for i = 1:numparts
    
    fprintf('Smoothing particle %d\n',i);
    tic
    particle = (psurf==i);
    [z,x,y] = ind2sub(size(particle),find(particle));
    %pprop = regionprops(particle,'Centroid');
    %center = [pprop.Centroid];
    
    %zcal corrects for the pixel compression from before
    y = y./zcal;
    x = x./zcal;
    z = (size(I,1)-z)./zcal;%1072 makes z-axis '0' the bottom of the box instead of the top
    center = mean([x y z]);
    
    if size(x,1)>200
        
    %map to spherical coordinates to perform smoothing
    xo = x-center(1);
    yo = y-center(2);
    zo = z-center(3);
    [phi,theta,rho]=cart2sph(xo,yo,zo);
    phi=phi(rho>0.8*median(rho));
    theta=theta(rho>0.8*median(rho));
    rho=rho(rho>0.8*median(rho));
    
    %remove outliers based on radius deviation from median
    dmedrho = rho-median(rho);
    phi=phi(dmedrho<0.3*max(dmedrho)&dmedrho>0.3*min(dmedrho));
    theta=theta(dmedrho<0.3*max(dmedrho)&dmedrho>0.3*min(dmedrho));
    rho=rho(dmedrho<0.3*max(dmedrho)&dmedrho>0.3*min(dmedrho));
    
    %generate periodic image in all directions
    ptheta = cat(1,theta-pi,theta-pi,theta-pi,theta,theta,theta,theta+pi,theta+pi,theta+pi);
    pphi = cat(1,phi-2*pi,phi,phi+2*pi,phi-2*pi,phi,phi+2*pi,phi-2*pi,phi,phi+2*pi);
    prho = cat(1,rho,rho,rho,rho,rho,rho,rho,rho,rho);
    
    %mesh and smooth
    p=-3*pi:pi/100:3*pi;t=-3*pi/2:pi/100:3*pi/2;
    [phiq,thetaq]=meshgrid(p,t);
    vq = griddata(pphi,ptheta,prho,phiq,thetaq,'cubic')+1.0/zcal;%plus 1.0 to account for watershed
    vq(isnan(vq))=nanmedian(nanmedian(vq));
    vq = smooth2a(vq,10,10);%requires smooth2a.m (by Greg Reeves)
    mesh(phiq,thetaq,vq);
    
    %removes the periodic images
    phiq = phiq(round(size(phiq,1)/3):round(2*size(phiq,1)/3)...
        ,round(size(phiq,2)/3):round(2*size(phiq,2)/3));
    thetaq = thetaq(round(size(thetaq,1)/3):round(2*size(thetaq,1)/3)...
        ,round(size(thetaq,2)/3):round(2*size(thetaq,2)/3));
    vq = vq(round(size(vq,1)/3):round(2*size(vq,1)/3)...
        ,round(size(vq,2)/3):round(2*size(vq,2)/3));
    
    %convert back to cartesian coordinates
    [xv,yv,zv]=sph2cart(phiq,thetaq,vq);
    toc
    %if you want to visualize particle mesh
     mesh(xv+center(1),yv+center(2),zv+center(3));axis equal

    totalpartvol = totalpartvol + parea(i);
    partscounted = partscounted + 1;

            if center(3)>size(I,2)/zcal/5 && center(3)<size(I,2)/zcal*4/5 && center(2)<size(II,3)/zcal*9/10
                %this removes image "particles" that are outside the box boundary 
                %(needs to be adjusted for your experiment)
                partcents(i,1)=center(1);
                partcents(i,2)=center(2);
                partcents(i,3)=center(3);

                partlocs(:,1,i)=reshape(phiq,numel(phiq),1);
                partlocs(:,2,i)=reshape(thetaq,numel(thetaq),1);
                partlocs(:,3,i)=reshape(vq,numel(vq),1);
            end
    end
end

partcents(any(partcents==0,2),:)=[];
partlocs(:,:,any(any(isnan(partlocs),1),2))=[];

numparts = size(partcents,1);

partdists = dist(partcents');
partdists(partdists>300) = 0;%identify possible contacts (<160 pixels apart) makes matrix sparse

contactphi = cell(numparts);
contacttheta = cell(numparts);
contactrad = cell(numparts);
%Identify contacts and calculate forces

parfor id = 1:numparts

    r1 = partlocs(:,3,id);
    phi1 = partlocs(:,1,id);
    tht1 = partlocs(:,2,id);
    fprintf('Looking for contacts on particle %d\n',id);
    
    for iid = 1:numparts
        
        distance = partdists(id,iid);
        
        if distance~=0
            
            r2 = partlocs(:,3,iid);
            
            if max(r1)+max(r2)>distance%possible contact?
                
                fprintf('Possible contact with particle %d\n',iid);
                
                phi2 = partlocs(:,1,iid)+pi;%invert particle 2
                tht2 = -partlocs(:,2,iid);
                
                %switch back to cartesian
                [x1,y1,z1] = sph2cart(phi1,tht1,r1);
                [x2,y2,z2] = sph2cart(phi2,tht2,r2);
                
                [phi1,tht1,r1] = cart2sph(x1,y1,z1);
                [phi2,tht2,r2] = cart2sph(x2,y2,z2);
                rtot = zeros(size(r1,1),1);
                rad2 = zeros(size(r1,1),1);
                
                for index = 1:size(x1,1)
                    for index2 = 1:size(x1,1)
                        if abs(phi2(index2)-phi1(index))<0.000001&&abs(tht2(index2)-tht1(index))<0.000001
                            rtot(index) = r1(index)+r2(index2);
                            rad2(index) = r2(index2);
                        end
                    end
                end
                
                dhat = (partcents(iid,:)-partcents(id,:))/distance;%d vector
                
                rad1 = r1;
                
                [xtot,ytot,ztot] = sph2cart(phi1,tht1,rtot);
                v = [xtot ytot ztot];
                
                theta = zeros(size(x1,1),1);
                dotprod = zeros(size(x1,1),1);
                
                for jj = 1:size(x1,1)
                    dotprod(jj) = sum(v(jj,:)/norm(v(jj,:)).*dhat(:)');
                    theta(jj) = acos(dotprod(jj));
                end
                
                v=v(theta<pi/2,:);
                rad1=rad1(theta<pi/2);
                rad2=rad2(theta<pi/2);

                dotprod = dotprod(theta<pi/2);
                
                [phi12,theta12,r12] = cart2sph(v(:,1),v(:,2),v(:,3));
                
                r12 = r12.*dotprod;
                rad1 = rad1.*dotprod;
                rad2 = rad2.*dotprod;
                phi12 = phi12(r12>distance);
                theta12 = theta12(r12>distance);
                
                contactphi{id,iid} = phi12;
                contacttheta{id,iid} = theta12;
                contactrad{id,iid} = [max(rad1) max(rad2)];
                if size(phi12,1)>0
                    fprintf('Particle %d contacted particle %d with angular range %d\n',id,iid,max(theta12)-min(theta12));
                end
            end
            
        end
    end
    
end

%Write particle and contact files for analysis
fprintf('Storing particle and contact files for scan %d\n',s);
tic
folder = sprintf('/filepath/particlefiles/scan%d',s);
if ~exist(folder, 'dir')
  mkdir(folder);
end
folder = sprintf('/filepath/contactfiles/scan%d',s);
if ~exist(folder, 'dir')
  mkdir(folder);
end

csvwrite(sprintf('/filepath/particlefiles/scan%d/centers%d.csv',s,s),partcents);

for i=1:numparts
    csvwrite(sprintf('/filepath/particlefiles/scan%d/particle%d.csv',s,i),partlocs(:,:,i));
end
file = sprintf('/filepath/contactfiles/scan%d/contactphi%d.mat',s,s);
save(file,'contactphi');
file = sprintf('/filepath/contactfiles/scan%d/contacttheta%d.mat',s,s);
save(file,'contacttheta');
file = sprintf('/filepath/contactfiles/scan%d/contactrad%d.mat',s,s);
save(file,'contactrad');
toc

file = sprintf('/filepath/contactfiles/scan%d/partdists%d.mat',s,s);
save(file,'partdists');
toc

%---------------------------------------
%begin calculating forces
%----------------------------------------

fprintf('Calculating forces and stress+fabric tensors for scan %d\n',s);
tic

pixeltometer = 0.243/2268;%estimated conversion factor

forces = zeros(numparts);
totalforce = zeros(numparts,1);
avgforce = zeros(numparts,1);
avgZ = zeros(numparts,1);
avgdelta = zeros(numparts,2);
avgrad = zeros(numparts,2);
centers = cell(numparts);
fabric = [0 0 0;0 0 0;0 0 0];
stress = [0 0 0;0 0 0;0 0 0];
branches = cell(numel(find(~cellfun(@isempty,contactphi))),1);
dex = 1;

contactarea = zeros(numparts);

for id1 = 1:numparts
    Z = 0;
    for id2 = 1:numparts
        if ~isempty(contactphi{id1,id2})
            
            phi = contactphi{id1,id2};
            theta = contacttheta{id1,id2};
            rad = contactrad{id1,id2};
            [xm,ym,zm] = sph2cart(phi,theta,sum(rad));
            d = partdists(id1,id2);
            
            avgrad(id1,id2) = sum(rad)/2;

            delta = sum(rad)-d;
            
            if rad(1) < 160 && rad(2) < 160 && delta > 0 && rad(1) > 80 && rad(2) > 80%radii cannot be too large
            
                avgdelta(id1,id2)=delta;
                [x1,y1,z1] = sph2cart(phi,theta,rad(1,1)-delta/2);%contact center
                center = mean([x1 y1 z1]);%contact center
                centers{id1,id2} = center;
                branch = center./norm(center);
                branches{dex,1} = branch;
                dex = dex+1;
                
                xyzdist = dist([x1,y1,z1]');

                %calculate contact force from Hertzian contact mechanics
                radcurve = 1/(1/rad(1)+1/rad(2));
                forces(id1,id2) = 4/3*1/(1.5/33750)*radcurve^0.5*delta^1.5*pixeltometer^2;

                %get fabric and stress tensors
                if forces(id1,id2) > 0 && forces(id1,id2) < 0.075 %force thresholds - adjust for your experiment
                    fabric = fabric + branch'*branch;
                    stress = stress + (branch'*branch).*forces(id1,id2)./(boxheight*0.18*0.243);%boxvolume*conversion factor
                    Z = Z+1;
                end
            end
        end
        if forces(id1,id2) > 0 && forces(id1,id2) < 0.075 %removing erroneous forces that are to low and too high (bad particles)
            totalforce(id1) = totalforce(id1)+forces(id1,id2);
        end
    end
    if Z~=0
        avgforce(id1) = totalforce(id1)/Z;
    end
    avgZ(id1) = Z;
end
toc
avgdelta(avgdelta==0) = NaN;

csvwrite(sprintf('/filepath/fabric/fabric%d.csv',s),fabric);
csvwrite(sprintf('/filepath/stress/stress%d.csv',s),stress);
csvwrite(sprintf('/filepath/forces/forces%d.csv',s),forces);
csvwrite(sprintf('/filepath/numparts/numparts%d.csv',s),numparts2);

end
