//iteratively spread out vert. & horizontal velocities among neighboring cells. We don't 
//average out over diagonal cells because horizontal velocities can only spread out to cells left or right
//and vertical velocities can only spread to cells up or down
 void lin_solve(int b, float[] x, float[] x0, float a, float c)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a*(    x[IX(i+1, j )]
                                    +x[IX(i-1, j  )]
                                    +x[IX(i  , j+1)]
                                    +x[IX(i  , j-1 )]
                                   
                           )) * cRecip;
                }
            }
    }
        //reset boundary cells to simulate reflection of particles off walls
        set_bnd(b, x);
    }
    
   //Function to diffuse the dye & velocity com ponents
  void diffuse(int b,float[] x,float[] x0,float diff,float dt)
  {
      //really dont have a clue why
      float a = dt * diff * (N-2) * (N-2);
      //spread out the diffusion over several iterations every render frame
      lin_solve(b,x,x0,a,1+6*a);
  
  }
  
 // try this after sim is set up : velocity of right & left walls for Vy stays same
 //                                velocity of up and down walls for Vy reflects
 //                                velocity of right & left walls for Vx reflects
 //                                velocity of up and down walls for Vx stays same
  void set_bnd(int b,float[] x)
  {
    
   //For dealing with vertical Vy Velocity and mirroring border cells
   
        for(int i = 1; i < N - 1; i++) {
            
            x[IX(i, 0 )] = b == 2 ? -x[IX(i, 1 )] : x[IX(i, 1 )];
            x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
        }
    
    //For dealing with horizontal Vx Velocity and mirroring border cells

        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
            x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
        }
    
   /* if(b == 2){
    for(int j=0;j<N-1;j++)
    {
      //set for Vy array
     
      x[IX(0, j ) ] =  -x[IX(1  , j)]; 
      x[IX(N-1, j)] =  -x[IX(N-2, j)];
        
    }
    }
    if(b == 1){
      for(int i=0; i<N-1;i++){
      x[IX(i,0)] =  -x[IX(i, 1 )];
      x[IX(i, N-1)] = -x[IX(i, N-2)];
      
      }
    }
    */
    //accounting for corner cells by taking average of its adjacent cells AFTER they are mirrored
    x[IX(0, 0)]       = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)]     = 0.5f * (x[IX(1, N-1)]+ x[IX(0, N-2)]);
    x[IX(N-1,0)]      = 0.5f * (x[IX(N-2,0)]+x[IX(N-1,1)]);
    x[IX(N-1,N-1)]    = 0.5f * (x[IX(N-1,N-2)] + x[IX(N-2,N-1)]);
    
  }
  
  
  //project method ensures the fluid velocity field remains incompressible
 void project(float[] velocX, float[] velocY, float [] p, float[] div)
{
    
  //this nested loop calculates divergence using the central difference method
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                         velocX[IX(i+1, j)]
                        -velocX[IX(i-1, j)]
                        +velocY[IX(i  , j+1)]
                        -velocY[IX(i  , j-1)]                        
                    )/N;        //Try with * to see if it produces a better effect
                p[IX(i, j)] = 0;//initialize pressure fields with 0
            }
        }
    
    set_bnd(0, div); 
    set_bnd(0, p);  //try commenting this and see effect
    lin_solve(0, p, div, 1, 6); // this solves for the pressure gradient
    
    
    //Subtract pressure field from velocity field to make velocity divergence free
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N;
            }
        }
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
    
}

 void advect(int b, float[] d, float[] d0,  float[] velocX, float[] velocY,float dt)
{
    float i0, i1, j0, j1;
    
    //scale time step by grid size
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
   
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    

        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
              
              // d = vt
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j )];
               
               
               //find the previous location of the field
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
   
                //ensure that x & y are within grid bounds to avoid "out of bounds" error
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floor(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floor(y);
                j1 = j0 + 1.0f; 
                
                
                s1 = x - i0; //extract the fractional difference for previous position (always < 1.0f)
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                
                int i0i = (int)i0; //x value of previous cell location
                int i1i = (int)i1; //x+1 value of previous cell location
                int j0i = (int)j0; //y value of previous cell location
                int j1i = (int)j1; //y+1 value of previous cell location
               
               
               /*bilinear interpolation to assign the velocity at X,Y to current velocity by considering the four neighbourhood 
               previous velocities  and their weights based on proximity of V(X,Y) to each grid point velocities */
                d[IX(i, j)] = 
                
                    s0 * (( t0 * d0[IX(i0i, j0i)])
                        +( t1 * d0[IX(i0i, j1i)]))
                                
                   +s1 * (( t0 * d0[IX(i1i, j0i)])
                        +( t1 * d0[IX(i1i, j1i)]));
            }
        }
    //set boundary velocities for current (and updated) velocities
    set_bnd(b, d);
}

void moveParticles(float[]Vx,float[]Vy,float[]ParticlesX,float[]ParticlesY,float dt)
{

  for(int i=0;i<N-1;i++)
  {
    for(int j=0;j<N-1;j++)
    {
      
      ParticlesX[IX(i,j)] += floor(Vx[IX(i,j)]*dt*pow(10,5));
      ParticlesY[IX(i,j)] += floor(Vy[IX(i,j)]*dt*pow(10,4));
     
     if(ParticlesX[IX(i,j)] < 0)
        ParticlesX[IX(i,j)] = 0;
     if(ParticlesX[IX(i,j)] > (N-1)*SCALE)
         ParticlesX[IX(i,j)] = N*SCALE;
     
     if(ParticlesY[IX(i,j)] < 0)
        ParticlesY[IX(i,j)] = 0;
     if(ParticlesY[IX(i,j)] > (N-1)*SCALE)
         ParticlesY[IX(i,j)] = N*SCALE;
    
    }
  
  }
  
  set_bnd(2,ParticlesY);
  set_bnd(1,ParticlesX);
  
  
}
