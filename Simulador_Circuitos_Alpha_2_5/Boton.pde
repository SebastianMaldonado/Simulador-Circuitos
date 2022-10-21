/*
|====================================================================|
*                               |Botón|         
*                            [Clase Padre]
* Descripción:                                                        
*   Clase definida para crear una interfaz donde poder ejecutar una
*   acción al ser presionado, se puede mostrar y devuelve un booleano
*   al pasar el ratón sobre él
|====================================================================|
*/
public class Boton {
  private int x, y;              //Posición del botón
  private int vx, vy;            //Posición virtual del botón (para movimiento del botón)
  private int tx, ty;            //Tamaño del botón
  private PImage imagen;         //Imagen del botón
  private PImage imagen_bloq;    //Imagen del botón bloqueado
  private PImage imagen_selc;    //Imagen del botón seleccionado
  
  
  public Boton (int x, int y, int tx, int ty) {
    this.x = x;
    this.y = y;
    this.vx = x;
    this.vy = y;
    this.tx = tx;
    this.ty = ty;
  }
  
  
  //-------------------------|Ingresar imagen|-------------------------\\
  public void ingresarImagen (PImage normal, PImage bloq, PImage selc) {
    this.imagen = normal;
    this.imagen_bloq = bloq;
    this.imagen_selc = selc;
    
    ajustarImagen(this.imagen);
    ajustarImagen(this.imagen_bloq);
    ajustarImagen(this.imagen_selc);
  }
  
  //
  private void ajustarImagen (PImage imagen) {
    if (imagen == null)
      return;
    
    imagen.resize(this.tx, this.ty);
  }
  
  //-------------------------|Sobre Caja|-------------------------\\
  public boolean sobreBoton (float ratonX, float ratonY) {
    if ((ratonX > this.vx) && (ratonX < this.vx + this.tx) && (ratonY > this.vy) && (ratonY < this.vy + this.ty))
      return true;    //Si el cursor está por encima de la ventana
    
    return false;   //Si el cursor NO está por encima de la ventana
  }
  
  //-------------------------|Mostrar Botón|-------------------------\\
  //Ingresar valores opciones para cargar una posición
  //No ingresar nada para ejecutarlo con la posición original del botón
  public void mostrar (PImage imagen, float x, float y) {
    if (imagen == null) {      //Si NO tiene una imagen
      rect (x, y, this.tx, this.ty);
      return;
    }                     
    
    //Si tiene una imagen
    image (imagen, x, y);
  }
  
  //-------------------------|Mostrar Botón Normal|-------------------------\\
  public void mostrar () {
    this.mostrar(this.imagen, this.x, this.y);
  }
  
  public void mostrar (float x, float y) {
    this.mostrar(this.imagen, x, y);
  }
  
  //-------------------------|Mostrar Botón Bloqueado|-------------------------\\
  public void mostrar_bloq () {
    fill(150);
    this.mostrar(this.imagen_bloq, this.x, this.y);
    fill(255);
  }
  
  public void mostrar_bloq (float x, float y) {
    fill(150);
    this.mostrar(this.imagen_bloq, x, y);
    fill(255);
  }
  
  //-------------------------|Mostrar Botón Seleccionado|-------------------------\\
  public void mostrar_selc () {
    fill(200);
    this.mostrar(this.imagen_selc, this.x, this.y);
    fill(255);
  }
  
  public void mostrar_selc (float x, float y) {
    fill(200);
    this.mostrar(this.imagen_selc, x, y);
    fill(255);
  }
}
