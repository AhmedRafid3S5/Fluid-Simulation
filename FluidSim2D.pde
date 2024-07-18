Fluid fluid;

//setting up window size etc
void settings(){
size(N*SCALE,N*SCALE);

}

//runs once at the start
void setup(){
 fluid = new Fluid(0.10,0, 0.0000001); //Fluid(timeStep,diffusion,viscosity)
}

void mouseDragged(){
  fluid.addDensity(mouseX/SCALE,mouseY/SCALE,100);
  float amtX = mouseX - pmouseX;
  float amtY = mouseY - pmouseY;
  fluid.addVelocity(mouseX/SCALE,mouseY/SCALE,amtX,amtY);
 
}


//continuously loops during the runtime
void draw() {
  background(0);
  fluid.renderD();
  //fluid.renderP();

  fluid.step(); // Perform a simulation step
}
