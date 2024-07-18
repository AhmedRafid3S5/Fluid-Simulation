//Declare Globals
final int N = 128;
final int iter = 10;
final int SCALE = 5;

//return the 1D equivalent index of a 2D array
int IX(int x, int y)
{
  x = constrain(x,0,N-1);
  y = constrain(y,0,N-1);
  return x + y*N; 
}


class Fluid{
    //Fluid cell properties
    int size;
    float dt; //the timeStep per simulation or Delta time
    float diff; //diffusion amount
    float visc; //Vescosity of the fluid
    
    float[] s; //previous density?
    float[] density; //density of each particle in the Fluid
    
    float[] Vx; //Horizontal velocity componenet in each cell
    float[] Vy; //Verticle velocity component in each cell

    float[] Vx0; //Previous Vx
    float[] Vy0; //Previous Vy
    
    float[] ParticlesX;
    float[] ParticlesY;
    
    
    
    //Constructor
    Fluid(float dt,float diffusion,float viscosity)
    {
      this.size = N;
      this.dt = dt;
      this.diff = diffusion;
      this.visc = viscosity;
      
      this.s = new float[N*N];
      this.density = new float[N*N];
      
      this.Vx = new float[N*N];
      this.Vy = new float[N*N];
      
      this.Vx0 = new float[N*N];
      this.Vy0 = new float[N*N];
      
      this.ParticlesX = new float[N*N];
      this.ParticlesY = new float[N*N];
      
     
    
    }
    
    void step()
{
    float visc     = this.visc;
    float diff     = this.diff;
    float dt       = this.dt;
    float[]Vx      = this.Vx;
    float[]Vy      = this.Vy;
    float[]Vx0     = this.Vx0;
    float[]Vy0     = this.Vy0;
    float[]s       = this.s;
    float[]density = this.density;
    float[]ParticlesX = this.ParticlesX;
    float[]ParticlesY = this.ParticlesY;
    
    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);
    
    
    project(Vx0, Vy0,Vx, Vy);
    moveParticles(Vx,Vy,ParticlesX,ParticlesY,dt);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);   
    moveParticles(Vx,Vy,ParticlesX,ParticlesY,dt);

    project(Vx, Vy, Vx0, Vy0);
    moveParticles(Vx,Vy,ParticlesX,ParticlesY,dt);
    
    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy,dt);
   
    
  println("Vx0[middle]: " + Vx0[IX(N*SCALE/2,N*SCALE/2)] + ", Vy0[0]: " + Vy0[IX(N*SCALE/2,N*SCALE/2)]);
}
    
    //a function to add some density to a particle/cell
    void addDensity(int x,int y, int amount)
    {
      int index = IX(x,y);
      this.density[index] += amount;
    }
    
    void addVelocity(int x,int y,float x_amount,float y_amount)
    {
      int index = IX(x,y);
      this.Vx[index] += x_amount;
      this.Vy[index] += y_amount;
      
    }
  
    void renderD(){
      
      colorMode(HSB,255);
      for(int i=0;i<N;i++){
         for(int j=0;j<N;j++)
         {
          float x = i * SCALE;
          float y = j * SCALE;
          float d = this.density[IX(i,j)];
        
          fill(d,255,255);
          noStroke();
          square(x,y,SCALE);
         }
      }
    
    
    }
    
    void renderP(){
    
    colorMode(HSB,255);
      for(int i=0;i<N;i++){
         for(int j=0;j<N;j++)
         {
          float x = i * SCALE;
          float y = j * SCALE;
          float d = this.density[IX(i,j)];
          float pX = this.ParticlesX[IX(i,j)];
          float pY = this.ParticlesY[IX(i,j)];
        
          fill(d,255,255);
          circle(pX,pY,SCALE);
         }
      }
    }
    
    void fadeDye(){
      for(int i=0;i<this.density.length;i++)
      {
         float d = density[i];
         density[i] = constrain(d-1,0,255);
      }
    }
    
    
    float[] velocityMagnitude() {
  float[] magnitude = new float[N*N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = IX(i, j);
      float vx = this.Vx[index];
      float vy = this.Vy[index];
      magnitude[index] = sqrt(vx * vx + vy * vy);
    }
  }
  return magnitude;
}

    
  
}
