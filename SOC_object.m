classdef SOC_object
    properties
        %parameters of the coil
        Rho = 1/0.28*1.68e-8 ;%ohm meter
        N = 100;%number of windings of coil 1
        H = 1e-3;%m height of coil 1
        D1 = 1e-3;%m inner diameter of coil 1
        D2 = 1.5e-3;%m out diameter of coil 1
        xpos = 1e-3;%distance between nozzle and coil in x
        zpos = 1e-3;%distance between nozzle and coil in z
        
        %"mesh" parameters
        nx = 11;%number of elements on the coil in the x direction
        nz = 10;%number of elements on the coil in the z direction
        nphi = 20;%number of elements of the coil in radial direction
        nxb = 20;%number of elements of the vertical boundaries
        nzb = 20;%number of elements of the horizonal boundaries
        is_coil = 0;%indicates the object is a coil (enables exta plotting options)
        is_circular = 0;%tell the solver to use circular elements, for comparing with Maxwell's equations
        
        ntot = 0;%total number of points on the mesh of the coil
        loops = 0;%loops of wire of the coil
        r = 0;%will contain the a vector containing the point on the coil
        dl = 0;%will contain a vectors of the direction of the current of the points in r
        el = 0;%will contain a normalised vector of the direciton of the current of the points in r
        dS = 0;%will contain a surface of each point
        dV = 0;%will contain the volume related to each point in rL
        
        %the content of these is the same as the content of the variables
        %above however these variables will require the content to be
        %placed on the gpu
        r_gpu = 0;
        dl_gpu = 0;
        el_gpu = 0;
        dS_gpu = 0;
        dV_gpu = 0;
        
        A_int = 0 %this will contain the field inside the wires
        
        
        taper_width = 0; %distance between the outer radius of the nozzle and the point where the taper starts
        taper_angle = 0; %the angle of the taper of the nozzle
        hole_offset = 0; %offset of the hole inside the nozzle
        rotation = 0; %rotation of the nozzle in degree
        
    end
    methods
        %set the geometry of the object to that of a hollow cylinder
        function obj = set_coil_geometry(obj, H,D1,D2,N,Rho,xpos,zpos)
            obj.H = H;
            obj.D1 = D1;
            obj.D2 = D2;
            obj.N = N;
            obj.Rho = Rho;
            obj.xpos = xpos;
            obj.zpos = zpos;
        end
        
        %set the geometry of the object to that of a typical 3D printer nozzle
        function obj = set_nozzle_geometry(obj, H,D1,D2,N,Rho,xpos,zpos,taper_width,taper_angle,hole_offset,rotation)
            obj.H = H;
            obj.D1 = D1;
            obj.D2 = D2;
            obj.N = N;
            obj.Rho = Rho;
            obj.xpos = xpos;
            obj.zpos = zpos;
            obj.taper_width = taper_width;
            obj.taper_angle = taper_angle;
            obj.hole_offset = hole_offset;
            obj.rotation = rotation;
        end
        
        %set mesh of the object to be uniformly spaced
        function obj = set_mesh(obj, nx,nz,nphi)
            obj.nx = nx;
            obj.nz = nz;
            obj.nphi = nphi;
        end
        
        %configure the mesh of the object to only be on the boundary. This
        %is usefull if due to the skin effect the current inside the nozzle
        %is negelectable
        function obj = set_boundary_mesh(obj, nx,nz,nphi,nxB,nzB)
            obj.nx = nx;
            obj.nz = nz;
            obj.nphi = nphi;
            obj.nxb = nxB;
            obj.nzb = nzB;
        end

        %this plots the geometry of the nozzle
        function plot_geometry(obj,C,prefix)
            if obj.r == 0
               error('error: cant plot geometry. Make sure the coil has been build.') 
            end
            scatter3(obj.r(:,1)*prefix(1),obj.r(:,2)*prefix(2),obj.r(:,3)*prefix(3),C)
            daspect([1 1 1])
        end

        %this plots the vectors forming the surface of an object
        function plot_surface_vectors(obj,prefix)
            quiver3(obj.r(:,1)*prefix(1),obj.r(:,2)*prefix(2),obj.r(:,3)*prefix(3),obj.dl(:,1)*prefix(1),obj.dl(:,2)*prefix(2),obj.dl(:,3)*prefix(3),0 )
            daspect([1 1 1])
        end
        
        %this plot the direction of the current inside a coil
        function plot_current_direction(obj)
            if ~obj.is_coil
                error('can only plot current direction if the object is a coil')
            end  
            quiver3(obj.r(:,1),obj.r(:,2),obj.r(:,3),obj.el(:,1),obj.el(:,2),obj.el(:,3))
            xlabel('x position (m)')
            ylabel('y position (m)')
            zlabel('z position (m)')
            daspect([1 1 1])
        end
        
        %this calculates the self inductance of a coil.
        function L = calculate_self_inductance(obj)
            if ~obj.is_coil
                error('can only calculate self inductance direction if the object is a coil')
            end 
            f = waitbar(0,'Calculating self inductance matrix');
            O = obj.H*((obj.D2-obj.D1)/2);
            M_CC = obj.inductance_matrix(obj.r);
            L = (obj.N^2)/(obj.nx*obj.nz*O)*ones(1,3)*diag(obj.dl.'*M_CC*obj.el);
            close(f)
        end
        
        %this calculates the mutual inductance between two coils.
        function L = calculate_mutual_inductance(obj,nozzle)
            if ((~obj.is_coil) || (~obj.is_coil))
                error('can only calculate mutual inductance if both objects are a coil')
            end 
            f = waitbar(0,'Calculating self inductance matrix');
            O = obj.H*((obj.D2-obj.D1)/2);
            M_CN = nozzle.inductance_matrix(obj.r);
            L = (obj.N^2)/(obj.nx*obj.nz*O)*ones(1,3)*diag(obj.dl.'*M_CN*nozzle.el);
            close(f)

        end
        
        %this calculates the change in inductance due to the presense of
        %another object
        function dL = calculate_inductance_change(obj,nozzle,frequency)
            w = frequency*2*pi;

            dZ = obj.calculate_impedance_change(nozzle,frequency);
            dL = dZ/1i/w;
        end
        
        %this calculates the change in impedance  
        function dZ = calculate_impedance_change(obj,nozzle,frequency)
            w = frequency*2*pi;
            
            f = waitbar(0,'Calculating coupling between coil and object');
            Mcn = nozzle.inductance_matrix_gpu(obj.r_gpu,f);
            waitbar(0,f,'Calculating objects self inductance matrix');
            Mnn = nozzle.inductance_matrix_gpu(nozzle.r_gpu,f);
            waitbar(0,f,'Calculating coupling between object and coil');
            Mnc = obj.inductance_matrix_gpu(nozzle.r_gpu,f);
            
            waitbar(0,f,'Calculating inductance change');
            if obj.is_circular
                O = obj.H*((obj.D2-obj.D1)/2)*pi/4;
            else
                O = obj.H*((obj.D2-obj.D1)/2);
            end
            part1 = -(obj.N^2*w^2)/(obj.nx*obj.nz*O*nozzle.Rho)*ones(1,3);
            part2 = obj.dl.'*Mcn;
            part3 = eye(nozzle.ntot)-1i*w/nozzle.Rho*Mnn;
            part4 = Mnc*obj.el;
            part5 = part3\part4;
            dZ = part1*diag(part2*part5);
            dZ = gather(dZ);
            close(f)
        end
        
        %this calculates the resistance of a coil
        function R = resistance(obj)
            if ~obj.is_coil
                error('can only calculate self inductance direction if the object is a coil')
            end
            if obj.is_circular
             A = obj.H*((obj.D2-obj.D1)/2)*pi/4/obj.N;
            else
             A = obj.H*((obj.D2-obj.D1)/2)/obj.N;
            end
            labs = vecnorm(obj.dl,2,2);
            L = sum(labs)/obj.nx/obj.nz*obj.N;
            R = obj.Rho*L/A;
        end
        
        %this plot the current density induced in the nozzle by the coil
        function plot_induced_current_density_2d(obj,nozzle,I,frequency,angle,prefix,offset)
            w = frequency*2*pi;
            
            f = waitbar(0,'Calculating coupling between coil and object');
            waitbar(0,f,'Calculating objects self inductance matrix');
            Mnn = nozzle.inductance_matrix_gpu(nozzle.r_gpu,f);
            waitbar(0,f,'Calculating coupling between object and coil');
            Mnc = obj.inductance_matrix_gpu(nozzle.r_gpu,f);
            
            O = obj.H*((obj.D2-obj.D1)/2);
            part1 = eye(nozzle.ntot)-1i*w/nozzle.Rho*Mnn;
            part2 = 1i*w*obj.N/nozzle.Rho/O*Mnc*obj.el*I;
                        
            JN = part1\part2;
            JN = double(gather(JN));
            JN_magn = sqrt(JN(:,1).^2+JN(:,2).^2+JN(:,3).^2);
            nozzle.plot_internal_scalar(abs(JN_magn),angle,'Induced current density (Am^{-2})',prefix,offset)
            close(f)
        end
        
        %calculate the skin depth for the object at a specific frequency
        function depth = skin_depth(obj,f)
            u0 = 4*pi*1e-7;
            depth = sqrt(2*obj.Rho/u0/2/pi/f);
        end
        
        %build the mesh for a cylindrical coil.
        function obj = build_coil(obj)
            f = waitbar(0,'Building coil');
            
            obj.is_coil = 1;
            obj.is_circular = 0;
            
            %calculate coil position vector
            obj.ntot = obj.nz*obj.nx*obj.nphi;
            obj.loops = obj.nz*obj.nx;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);
            obj.dS = zeros(obj.ntot,1);

            dx = (obj.D2-obj.D1)/2/obj.nx;
            x = obj.D1/2+dx/2:dx:obj.D2/2;

            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;

            dz = obj.H/obj.nz;
            z = dz/2:dz:obj.H-dz/2;

            for i1 = 1:obj.nz
                waitbar(i1/obj.nz,f)
                for i2 = 1:obj.nx
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = x(i2)*cos(phi(i3))+obj.xpos;
                        obj.r(loc,2) = x(i2)*sin(phi(i3));
                        obj.r(loc,3) = z(i1)+obj.zpos;
                        
                        obj.dl(loc,1) = -x(i2)*sin(phi(i3))*dphi;
                        obj.dl(loc,2) = x(i2)*cos(phi(i3))*dphi;
                        obj.dl(loc,3) = 0;
                        
                        obj.el(loc,1) = -sin(phi(i3));
                        obj.el(loc,2) = cos(phi(i3));
                        obj.el(loc,3) = 0;
                        
                        obj.dS(loc) = ((x(i2)+dx/2)^2-(x(i2)-dx/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dz*obj.dS(loc);
                    end
                end
            end
            obj.r_gpu = gpuArray(single(obj.r));
            obj.dl_gpu = gpuArray(single(obj.dl));
            obj.el_gpu = gpuArray(single(obj.el));
            obj.dS_gpu = gpuArray(single(obj.dS));
            obj.dV_gpu = gpuArray(single(obj.dV));
            
            close(f)
            
        end
        
        %build the mesh for a single loop of circular wire
        function obj = build_single_loop(obj)
            f = waitbar(0,'Building coil');
            
            obj.is_coil = 1;
            obj.is_circular = 1;
            
            %calculate coil position vector
            obj.ntot = obj.nz*obj.nx*obj.nphi;
            obj.loops = obj.nz*obj.nx;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);
            obj.dS = zeros(obj.ntot,1);

            r_circ = (obj.D2-obj.D1)/4;
            r_avg = (obj.D2+obj.D1)/4;
            

            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;

            dz = obj.H/obj.nz;
            z = dz/2:dz:obj.H-dz/2;

            for i1 = 1:obj.nz
                waitbar(i1/obj.nz,f)
                for i2 = 1:obj.nx
                    
                    y_circ = abs(obj.H/2-z(i1));
                    x_circ = sqrt(r_circ^2-y_circ^2);
                    dx = 2*x_circ/obj.nx;
                    x = r_avg-x_circ+dx/2:dx:r_avg+x_circ-dx/2;
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = x(i2)*cos(phi(i3))+obj.xpos;
                        obj.r(loc,2) = x(i2)*sin(phi(i3));
                        obj.r(loc,3) = z(i1)+obj.zpos;
                        
                        obj.dl(loc,1) = -x(i2)*sin(phi(i3))*dphi;
                        obj.dl(loc,2) = x(i2)*cos(phi(i3))*dphi;
                        obj.dl(loc,3) = 0;
                        
                        obj.el(loc,1) = -sin(phi(i3));
                        obj.el(loc,2) = cos(phi(i3));
                        obj.el(loc,3) = 0;
                        
                        obj.dS(loc) = ((x(i2)+dx/2)^2-(x(i2)-dx/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dz*obj.dS(loc);
                    end
                end
            end
            obj.r_gpu = gpuArray(single(obj.r));
            obj.dl_gpu = gpuArray(single(obj.dl));
            obj.el_gpu = gpuArray(single(obj.el));
            obj.dS_gpu = gpuArray(single(obj.dS));
            obj.dV_gpu = gpuArray(single(obj.dV));
            
            close(f)
            
        end
        
        %build build a mesh that only has points at the boundary of hollow
        %cylinder
        function obj = build_coil_surface(obj,depth)
            f = waitbar(0,'Building surface coil');
            
            obj.is_coil = 0;
            obj.is_circular = 0;
            
            %calculate coil position vector
            obj.ntot = obj.nzb*obj.nx*obj.nphi+obj.nz*obj.nxb*obj.nphi;
            obj.loops = 0;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);
            obj.dS = zeros(obj.ntot,1);

            dxb = 2*depth/obj.nxb;
            n1 = floor(obj.nxb/2);
            n2 = ceil(obj.nxb/2);
            x1 = obj.D1/2:dxb:obj.D1/2+dxb*(n1-1);
            x2 = obj.D2/2-dxb*(n2-1):dxb:obj.D2/2;
            xb = [x1,x2];
            
            dx = (obj.D2-obj.D1)/2/obj.nx;
            x = obj.D1/2+dx/2:dx:obj.D2/2-dx/2;
            
            
            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;
            
            dzb = 2*depth/obj.nzb;
            n1 = floor(obj.nzb/2);
            n2 = ceil(obj.nzb/2);
            z1 = 0:dzb:dzb*(n1-1);
            z2 = obj.H-dzb*(n2-1):dzb:obj.H;
            zb = [z1,z2];
            
            dz = (obj.H-2*depth)/(obj.nz+1);
            z = depth+dz:dz:obj.H-depth-dz;

            for i1 = 1:obj.nxb
                waitbar(i1/obj.nxb/2,f)
                for i2 = 1:obj.nz
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nz*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = xb(i1)*cos(phi(i3))+obj.xpos;
                        obj.r(loc,2) = xb(i1)*sin(phi(i3));
                        obj.r(loc,3) = z(i2)+obj.zpos;
                        
                        obj.dl(loc,1) = -xb(i1)*sin(phi(i3))*dphi;
                        obj.dl(loc,2) = xb(i1)*cos(phi(i3))*dphi;
                        obj.dl(loc,3) = 0;
                        
                        obj.el(loc,1) = -sin(phi(i3));
                        obj.el(loc,2) = cos(phi(i3));
                        obj.el(loc,3) = 0;
                        
                        obj.dS(loc) = ((xb(i1)+dxb/2)^2-(xb(i1)-dxb/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dz*obj.dS(loc);
                    end
                end
            end
            
            for i1 = 1:obj.nzb
                waitbar(i1/obj.nzb/2+0.5,f)
                for i2 = 1:obj.nx
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3+obj.nz*obj.nxb*obj.nphi;
                        obj.r(loc,1) = x(i2)*cos(phi(i3))+obj.xpos;
                        obj.r(loc,2) = x(i2)*sin(phi(i3));
                        obj.r(loc,3) = zb(i1)+obj.zpos;
                        
                        obj.dl(loc,1) = -x(i2)*sin(phi(i3))*dphi;
                        obj.dl(loc,2) = x(i2)*cos(phi(i3))*dphi;
                        obj.dl(loc,3) = 0;
                        
                        obj.el(loc,1) = -sin(phi(i3));
                        obj.el(loc,2) = cos(phi(i3));
                        obj.el(loc,3) = 0;
                        
                        obj.dS(loc) = ((x(i2)+dx/2)^2-(x(i2)-dx/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dzb*obj.dS(loc);
                    end
                end
            end
            
            obj.r_gpu = gpuArray(single(obj.r));
            obj.dl_gpu = gpuArray(single(obj.dl));
            obj.el_gpu = gpuArray(single(obj.el));
            obj.dS_gpu = gpuArray(single(obj.dS));
            obj.dV_gpu = gpuArray(single(obj.dV));
            close(f)
            
        end
        
        %build a mesh with only has points at the boundary of a 3d printer
        %nozzle shaped object
        function obj = build_nozzle_surface(obj,depth)
            f = waitbar(0,'Building surface coil');
            
            obj.is_coil = 0;
            obj.is_circular = 0;
            %calculate coil position vector
            obj.ntot = obj.nzb*obj.nx*obj.nphi+obj.nz*obj.nxb*obj.nphi;
            obj.loops = 0;
            obj.r = zeros(obj.ntot,3);
            obj.dl = zeros(obj.ntot,3);
            obj.dV = zeros(obj.ntot,1);
            obj.dS = zeros(obj.ntot,1);
            
            dphi = 2*pi/obj.nphi;
            phi = dphi/2:dphi:2*pi;
            
            dzb = 2*depth/obj.nzb;
            n1 = floor(obj.nzb/2);
            n2 = ceil(obj.nzb/2);
            z1 = 0:dzb:dzb*(n1-1);
            z2 = obj.H-dzb*(n2-1):dzb:obj.H;
            zb = [z1,z2];
            
            dz = (obj.H-2*depth)/obj.nz;
            z = depth+dz/2:dz:obj.H-depth-dz/2;
               
            for i3 = 1:obj.nphi
                if phi(i3) <= pi  
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(phi(i3)))^2);
                else
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(2*pi-phi(i3)))^2);
                end
            end
            
            %calculate x locations for the mesh on the left and right side of the
            %nozzle mesh
            for i2 = 1:obj.nz
                for i3 = 1:obj.nphi
                    taper_height = obj.taper_width*tan(obj.taper_angle);
                    
                    dxb = 2*depth/obj.nxb;
                    n1 = floor(obj.nxb/2);
                    n2 = ceil(obj.nxb/2);
                    
                    if z(i2) > taper_height
                        x1 = xstart(i3):dxb:xstart(i3)+dxb*(n1-1);
                        x2 = obj.D2/2-dxb*(n2-1):dxb:obj.D2/2;
                        xb{i2,i3} = [x1,x2];

                    else
                        xstop = obj.D2/2-obj.taper_width/taper_height*(taper_height-z(i2));
                        
                        x1 = xstart(i3):dxb:xstart(i3)+dxb*(n1-1);
                        x2 = xstop-dxb*(n2-1):dxb:xstop;
                        xb{i2,i3} = [x1,x2];
                    end
                end
            end
            
            %calculate r,dl,el,ds,dV for the left and right side of the
            %mesh
            for i1 = 1:obj.nxb
                waitbar(i1/obj.nxb/2,f)
                for i2 = 1:obj.nz
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nz*obj.nphi+(i2-1)*obj.nphi+i3;
                        obj.r(loc,1) = xb{i2,i3}(i1)*cos(phi(i3));
                        obj.r(loc,2) = xb{i2,i3}(i1)*sin(phi(i3));
                        obj.r(loc,3) = z(i2);
                        
                        obj.dS(loc) = ((xb{i2,i3}(i1)+dxb/2)^2-(xb{i2,i3}(i1)-dxb/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dz*obj.dS(loc);
                    end
                end
            end
            
            %calculate x locations for the mesh on the top and bottom side of the
            %nozzle mesh
            for i3 = 1:obj.nphi
                if phi(i3) <= pi  
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(phi(i3)))^2);
                else
                    xstart(i3) = obj.hole_offset*sin(phi(i3))+sqrt((obj.D1/2)^2 - (obj.hole_offset*cos(2*pi-phi(i3)))^2);
                end
            end
            
            for i1 = 1:obj.nzb
                for i3 = 1:obj.nphi
                    taper_height = obj.taper_width*tan(obj.taper_angle);
                    
                    
                    dxb = 2*depth/obj.nxb;
                    n1 = floor(obj.nxb/2);
                    n2 = ceil(obj.nxb/2);
                    
                    if zb(i1) >= taper_height

                        dx{i1,i3} = (obj.D2/2-xstart(i3))/obj.nx;
                        x{i1,i3} = xstart(i3)+dx{i1,i3}/2:dx{i1,i3}:obj.D2/2-dx{i1,i3}/2;
                    else
                        xstop = obj.D2/2-obj.taper_width/taper_height*(taper_height-zb(i1));
                        

                        dx{i1,i3} = (xstop-xstart(i3))/obj.nx;
                        x{i1,i3} = xstart(i3)+dx{i1,i3}/2:dx{i1,i3}:xstop;
                    end
                end
            end
            
            %calculate r,dl,el,ds,dV for the top and bottom side of the
            %mesh
            for i1 = 1:obj.nzb
                waitbar(i1/obj.nzb/2+0.5,f)
                for i2 = 1:obj.nx
                    for i3 = 1:obj.nphi
                        loc = (i1-1)*obj.nx*obj.nphi+(i2-1)*obj.nphi+i3+obj.nz*obj.nxb*obj.nphi;
                        obj.r(loc,1) = x{i1,i3}(i2)*cos(phi(i3));
                        obj.r(loc,2) = x{i1,i3}(i2)*sin(phi(i3));
                        obj.r(loc,3) = zb(i1);
                        
                        obj.dS(loc) = ((x{i1,i3}(i2)+dx{i1,i3}/2)^2-(x{i1,i3}(i2)-dx{i1,i3}/2)^2)*pi/obj.nphi;
                        obj.dV(loc) = dzb*obj.dS(loc);
                    end
                end
            end
            
            

            rtemp = obj.r(:,1);
            obj.r(:,1) = obj.r(:,2);
            obj.r(:,2) = rtemp;
            
            rotation_rad =  obj.rotation*pi/180;
            obj.r(:,1) = cos(rotation_rad)*obj.r(:,1)-sin(rotation_rad)*obj.r(:,3);
            obj.r(:,3) = sin(rotation_rad)*obj.r(:,1)+cos(rotation_rad)*obj.r(:,3);
            
            obj.r(:,1) = obj.r(:,1) + obj.xpos;
            obj.r(:,3) = obj.r(:,3) + obj.zpos;
            
            obj.r_gpu = gpuArray(single(obj.r));
            obj.dS_gpu = gpuArray(single(obj.dS));
            obj.dV_gpu = gpuArray(single(obj.dV));
            close(f)
            
        end
    end
    methods (Access = private)
        function L = inductance_matrix(obj,rd)
            u0 = 4*pi*1e-7;
            nd = length(rd);
            L = zeros(nd,obj.ntot);
            parfor i1 = 1:nd
                rd_m = ones(obj.ntot,1)*rd(i1,:);
                r_dif = rd_m-obj.r;
                r_abs = vecnorm(r_dif,2,2);   
                L(i1,:) = u0/4/pi./r_abs.*obj.dV;
            end
            L(isinf(L)) = 0;
        end
        
        function L = inductance_matrix_gpu(obj,rd,f)
            u0 = gpuArray(single(4*pi*1e-7));
            nd = length(rd);
            L = gpuArray(single(zeros(nd,obj.ntot)));
            for i1 = 1:nd
                if mod(i1,100) == 0
                    waitbar(i1/nd,f)
                end
                rd_m = gpuArray(single(ones(obj.ntot,1)))*rd(i1,:);
                r_dif = rd_m-obj.r_gpu;
                r_abs = vecnorm(r_dif,2,2);   
                L(i1,:) = u0/4/pi./r_abs.*obj.dV_gpu;
            end
            L(isinf(L)) = 0;
        end
        
        function plot_trisurf(~,x,y,z)
            dt = delaunayTriangulation(x,y) ;
            tri = dt.ConnectivityList ;
            xi = dt.Points(:,1) ; 
            yi = dt.Points(:,2) ; 
            F = scatteredInterpolant(x,y,z);
            zi = F(xi,yi) ;
            trisurf(tri,xi,yi,zi) 
            view(2)
            shading interp
        end     
        
            
        
        function plot_internal_scalar(obj,phi,direction,label,prefix,offset)
            n = round(((unwrap(direction)+pi)/2/pi)*(obj.nphi-1))+1;
            use = n:obj.nphi:obj.ntot;
            magn = sqrt(obj.r(use,1).^2+obj.r(use,2).^2);
            z = obj.r(use,3);
            obj.plot_trisurf(magn*prefix(1)+offset(1),z*prefix(2)+offset(2),real(phi(use))*prefix(3)+offset(3));
            c = colorbar;
            c.Label.String = label;
        end
    end
end