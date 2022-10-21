/*
|====================================================================|
*                               |Nodo|         
* Descripción:                                                        
*   Clase utilizada para almacenar todos los ramales que se unen
*   en una intersección o nodo, contiene la posición, un nombre
*   y la lista de ramales que se conectan en él
|====================================================================|
*/

public class Nodo {
  private String nombre;              //Nombre del Nodo
  private int pos_x, pos_y;           //Posición del Nodo en el plano
  private ArrayList<Ramal> ramales;   //Ramales que conecta el Nodo
  
  
  public Nodo (int x, int y) {
    this.ramales = new ArrayList<Ramal>();
    this.nombre = "-";
    this.pos_x = x;
    this.pos_y = y;
  }
  
  //-------------------------|Agregar Ramal|-------------------------\\
  public void nombrar (String nombre) {
    this.nombre = nombre;
  }
  
  //-------------------------|Agregar Ramal|-------------------------\\
  public void agregarRamal (Ramal ramal) {
    this.ramales.add(ramal);
  }
  
  //-------------------------|Eliminar Ramal|-------------------------\\
  public void eliminarRamal (Ramal ramal) {
    for (int i = 0; i < this.ramales.size(); i++)    //Eliminar Ramal a Eliminar
      if (ramal == this.ramales.get(i)) {
        this.ramales.remove(i);
        return;
      }
  }
  
  //-------------------------|Mostrar Nodo|-------------------------\\
  public void mostrar (Motor_grafico motor) {
    motor.mostrar("Nodo", this.pos_x - 50, this.pos_y - 50, 100, 100);
  }
  
  //-------------------------|Mostrar Nodo Seleccionado|-------------------------\\
  public void mostrarAreaSeleccion (int rad, Motor_grafico motor) {
    motor.mostrarCirculo(this.pos_x, this.pos_y, rad);
  }
  
  //-------------------------|Siguiente Nodo|-------------------------\\
  public Nodo siguiente (Ramal ramal) {
    return ramal.siguienteNodo(this);
  }
  
  //-------------------------|Obtener Ramales|-------------------------\\
  public ArrayList<Ramal> obtenerRamales () {
    return this.ramales;
  }
  
  //-------------------------|Obtener Vecinos|-------------------------\\
  public ArrayList<Nodo> obtenerVecinos () {
    ArrayList<Nodo> vecinos = new ArrayList<Nodo>();
    
    for (int i = 0; i < this.ramales.size(); i++)
      vecinos.add(this.ramales.get(i).siguienteNodo(this));

    return vecinos;
  }
  
  //-------------------------|Obtener Posición en X|-------------------------\\
  public int obtenerPosX () {
    return this.pos_x;
  }
  
  //-------------------------|Obtener Posición en Y|-------------------------\\
  public int obtenerPosY () {
    return this.pos_y;
  }
  
  //-------------------------|Obtener Nombre|-------------------------\\
  public String obtenerNombre () {
    return this.nombre;
  }
  
  //-------------------------|Modificar Posiciones|-------------------------\\
  public void modificarPos (int x, int y) {
    int ant_x = this.pos_x;
    int ant_y = this.pos_y;
    
    this.pos_x = x;
    this.pos_y = y;
    
    for (int i = 0; i < this.ramales.size(); i++) {
      //Actualizar longitud
      this.ramales.get(i).calcularLongitud();
      
      //Si la longitud actual es menor que la mínima no modificar más
      //println(this.ramales.size());
      
      if (this.ramales.get(i).longitud() <= this.ramales.get(i).longitudMin())
        modificarPos(ant_x, ant_y);
    }
  }
}
