/*
|====================================================================|
*                              |Ramal|         
* Descripción:                                                        
*   Clase utilizada para almacenar toda la información y las 
*   funciones necesarias para el funcionamiento de cada ramal
|====================================================================|
*/

public class Ramal {
  
  //Composición del Ramal
    private Nodo nodo1, nodo2;                    //Nodos a los que conecta en cada extremo
    private ArrayList<Componente> componentes;    //Lista de Componentes con los que conecta
  
  //Propiedades del Ramal
    private float longitud;                       //Longitud del Ramal
    private float long_min;                       //Longitud Mínima a la que se puede reducir el Ramal
  
  //Variables físicas
    private float corriente;                      //Corriente que atravisa el Ramal
    private float voltaje;                        //Diferencia de Potencial del Ramal
    private float potencia;                       //Potencia eléctrica del Ramal
    private boolean orientacion;                  //Orientación de la corriente en el ramal
  
  
  
  public Ramal (Nodo nodo1, Nodo nodo2) {
    //Ingresar Nodos extremos
    this.nodo1 = nodo1;
    this.nodo2 = nodo2;
    
    calcularLongitud();
    this.long_min = 0;
    
    //Conectar Nodos mediante Ramal
    this.nodo1.agregarRamal (this);
    this.nodo2.agregarRamal (this);
    
    this.componentes = new ArrayList<Componente>();
  }
  
  
  //---------------------------------------------------------------------------\\
  //                          |Propiedades del Ramal|                          \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Obtener Componentes|-------------------------\\
  public ArrayList<Componente> obtenerComponentes () {
    return this.componentes;
  }
  
  //-------------------------|Eliminar Componente|-------------------------\\
  public void eliminarComponente (Componente componente) {
    for (int i = 0; i < this.componentes.size(); i++) {
      if (this.componentes.get(i) == componente) {
        this.componentes.remove(i);
        actualizarLongMin();
        return;
      }
    }
  }
  
  //-------------------------|Modificar Corriente|-------------------------\\
  public void modificarCorriente (float corriente) {
    if (corriente == 0) {
      this.corriente = 0;
      return;
    }
    
    this.corriente = this.corriente + corriente;
  }
  
  //-------------------------|Ajustar Dirección de la Corriente|-------------------------\\
  public void ajustarDirCorriente () {
    if (this.corriente < 0)
      this.orientacion = false;
    
    this.corriente = Calculadora_geo.valorAbs(this.corriente);
  }
  
  //-------------------------|Obtener Información|-------------------------\\
  public String[] obtenerInfo () {
    String[] info;
    if (this.obtenerResistencia() > 0)
      info = new String[] {this.obtenerNombre(), this.obtenerResistencia() + "", this.corriente + "", this.voltaje + "", this.potencia + "", "", ""};
    else
      info = new String[] {this.obtenerNombre(), "--", this.corriente + "", "--", "--", "", ""};
    
    if (this.obtenerFEM() > 0) {
      info[5] = this.obtenerFEM() + "";
      info[6] = (this.corriente * this.obtenerFEM()) + "";
      return info;
    }
    
    info[5] = "--";
    info[6] = "--";
    
    return info; 
  }
  
  //-------------------------|Obtener Resistencia|-------------------------\\
  public float obtenerResistencia () {
    return obtenerValorComponente("Resistor");
  }
  
  //-------------------------|Obtener Resistencia|-------------------------\\
  public float obtenerFEM () {
    return obtenerValorComponente("Fuente");
  }
  
  //-------------------------|Obtener Valor de Componentes|-------------------------\\
  private float obtenerValorComponente (String componente) {
    float res = 0;
    
    for (int i = 0; i < this.componentes.size(); i++) {
      if (this.componentes.get(i).obtenerNombre().equals(componente)) 
        res = this.componentes.get(i).arrojarValor() + res;
    }
    
    return res;
  }
  
  //-------------------------|Obtener Nodo 1|-------------------------\\
  public Nodo obtenerNodo1 () {
    return this.nodo1;
  }
  
  //-------------------------|Obtener Nodo 2|-------------------------\\
  public Nodo obtenerNodo2 () {
    return this.nodo2;
  }
  
  //-------------------------|Obtener Nombre|-------------------------\\
  public String obtenerNombre () {
    return this.nodo1.obtenerNombre() + this.nodo2.obtenerNombre();
  }

  //-------------------------|Obtener Longitude|-------------------------\\
  public float longitud () {
    return this.longitud;
  }
  
  //-------------------------|Obtener Nombre|-------------------------\\
  public float longitudMin () {
    return this.long_min;
  }
  
  //-------------------------|Modificar Posiciones|-------------------------\\
  public Nodo siguienteNodo (Nodo nodo) {
    if (nodo == nodo1)
      return nodo2;
      
    if (nodo == nodo2)
      return nodo1;
      
    return null;
  }

  //-------------------------|Agregar Componentes|-------------------------\\
  public void agregarComponente (Componente componente) {
    if (this.long_min < componente.obtenerPos() + 300)
      this.long_min = componente.obtenerPos() + 300;
    
    this.componentes.add(componente);
  }
  
  //-------------------------|Actualizar Longitud Minima|-------------------------\\
  public void actualizarLongMin () {
    float may = 0;
    for (int i = 0; i < this.componentes.size(); i++)
      if (may < this.componentes.get(i).obtenerPos() + 300)
        may = this.componentes.get(i).obtenerPos() + 300;
        
    this.long_min = may;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                          |Modificación del Ramal|                         \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Separar Ramal|-------------------------\\
  public void separarRamal (float pos, Nodo nuevo, Circuito circuito) {
    //Crear nuevos ramales
    Ramal ramal1 = new Ramal (nodo1, nuevo);
    Ramal ramal2 = new Ramal (nuevo, nodo2);
    
    //Separar componentes del ramal dividido
    for (int i = 0; i < this.componentes.size(); i++) {
      if (this.componentes.get(i).obtenerPos() < pos) {
        ramal1.agregarComponente(this.componentes.get(i));
      } else {
        this.componentes.get(i).modificarPos(this.componentes.get(i).obtenerPos() - pos);
        ramal2.agregarComponente(this.componentes.get(i));
      }
    }
    
    //Agregar Ramales nuevos al circuito y eliminar el antiguo
    circuito.agregarRamal(ramal1);
    circuito.agregarRamal(ramal2);
    circuito.eliminarRamal(this);
  }
  
  
  //---------------------------------------------------------------------------\\
  //                          |Funciones Geométricas|                          \\
  //---------------------------------------------------------------------------\\  
  
  //-------------------------|Calcular Longitud del Ramal|-------------------------\\
  private void calcularLongitud () {
    int x1 = this.nodo1.obtenerPosX();
    int y1 = this.nodo1.obtenerPosY();
    int x2 = this.nodo2.obtenerPosX();
    int y2 = this.nodo2.obtenerPosY();
    this.longitud = Calculadora_geo.calcularDistancia (x1, y1, x2, y2);
  }
  
  //-------------------------|Calcular Punto más cercano al Ramal|-------------------------\\
  public float[] puntoCercanoRamal (float px, float py) {
    float x1 =  this.nodo1.obtenerPosX();
    float y1 =  this.nodo1.obtenerPosY();
    float x2 =  this.nodo2.obtenerPosX();
    float y2 =  this.nodo2.obtenerPosY();
      
    return Calculadora_geo.puntoPerpendicular(px, py, x1, y1, x2, y2);
  }
  
  //-------------------------|Posición Componente en el Ramal|-------------------------\\
  public float calcularPosComponente (float px, float py) {
    float[] punto = puntoCercanoRamal(px, py);
    return Calculadora_geo.calcularDistancia(punto[0], punto[1], this.nodo1.obtenerPosX(), this.nodo1.obtenerPosY());
  }
  
  
  //---------------------------------------------------------------------------\\
  //                           |Utilidades Gráficas|                           \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Mostrar Componentes|-------------------------\\
  public void mostrar (Motor_grafico motor) {
    motor.mostrarLinea(nodo1.obtenerPosX(), nodo1.obtenerPosY(), nodo2.obtenerPosX(), nodo2.obtenerPosY());
    
    for (int i = 0; i < this.componentes.size(); i++) {
      int[] x = new int[] {this.nodo1.obtenerPosX(), this.nodo2.obtenerPosX()};
      int[] y = new int[] {this.nodo1.obtenerPosY(), this.nodo2.obtenerPosY()};
      this.componentes.get(i).mostrar(x, y, motor);
    }
  }
  
  //-------------------------|Mostrar Área de Selección|-------------------------\\
  public void mostrarAreaSeleccion (int rad, Motor_grafico motor) {
    int x1 = nodo1.obtenerPosX();
    int y1 = nodo1.obtenerPosY();
    int x2 = nodo2.obtenerPosX();
    int y2 = nodo2.obtenerPosY();
  
    int[][] area_cuadrada = Calculadora_geo.calcularAreaCuadrada(x1, y1, x2, y2, rad);
    
    motor.mostrarLinea(area_cuadrada[0][0], area_cuadrada[1][0], area_cuadrada[0][1], area_cuadrada[1][1]);
    motor.mostrarLinea(area_cuadrada[0][1], area_cuadrada[1][1], area_cuadrada[0][2], area_cuadrada[1][2]);
    motor.mostrarLinea(area_cuadrada[0][2], area_cuadrada[1][2], area_cuadrada[0][3], area_cuadrada[1][3]);
    motor.mostrarLinea(area_cuadrada[0][3], area_cuadrada[1][3], area_cuadrada[0][0], area_cuadrada[1][0]);
  }
}
