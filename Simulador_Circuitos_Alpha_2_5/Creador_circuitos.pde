/*
|====================================================================|
*                        |Creador de Circuitos|         
* Descripción:                                                        
*   Clase utilizada para la construcción del circuito según los
*   eventos que vaya ingresando el usuario por medio de la interfaz
|====================================================================|
*/

public class Creador_circuitos {
  
  //Modos de actividad
    private String modo;          //Modo de creación de componentes
  
  //Creación de cables
    private int cont;             //Contador de pasos para ingreso de componentes
    private int ult_x, ult_y;     //Coordenadas del último click registrado
    private boolean modificar;    //Booleano que permite mantener la modificación de un objeto
  
  //Creación de Componentes
    private int valor;            //Valor de la característica del componente
    private boolean orientacion;  //Orientación del componente
  
  //Últimos elementos de Selección
    private Nodo ult_nodo;        //Último nodo seleccionado
    private Ramal ult_ramal;      //Último ramal seleccionado
    private Componente ult_comp;  //Último componente seleccionado
  
  
  
  public Creador_circuitos () {
    this.modo = "Ninguno";
    this.cont = 0;
    this.ult_x = 0;
    this.ult_y = 0;
    this.valor = 0;
    this.orientacion = true;
    this.modificar = false;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                         |Actividad según Eventos|                         \\
  //---------------------------------------------------------------------------\\  
  
  //---------------------------|Modificar|---------------------------\\
  public void modificar (float ratonX, float ratonY, int rad, Selector selector) {
    if (!this.modificar)
      return;
      
    if (this.ult_nodo == null && this.ult_comp == null) {
      this.modificar = false;
      return;
    }
      
    if (!this.modo.equals("Seleccion")) {
      this.modificar = false;
      this.ult_nodo = null;
      return;
    }
      
    if (this.ult_nodo != null) {
      this.ult_nodo.modificarPos(Math.round(ratonX), Math.round(ratonY));
      return;
    }
    
    if (this.ult_ramal == null)
      return;
    
    //Calcular Nueva Posición del Componente
    int pos = Math.round(this.ult_ramal.calcularPosComponente(ratonX, ratonY));
    
    this.ult_ramal.eliminarComponente(this.ult_comp);
    
    //Validar distancia de modificación
    float[] punto =  this.ult_ramal.puntoCercanoRamal(ratonX, ratonY);
    float dist = Calculadora_geo.calcularDistancia(punto[0], punto[1], ratonX, ratonY);
    if (dist > rad)
      return;
    
    if (validarColocacionComponente (ratonX, ratonY, rad, selector, this.ult_ramal))
      this.ult_comp.modificarPos(pos);
    
    this.ult_ramal.agregarComponente(this.ult_comp);
  }
  
  //---------------------------|Accionar|---------------------------\\
  public void accionar (float ratonX, float ratonY, int rad, Circuito circuito, Selector selector) {    
    if (this.modo.equals("Ninguno"))        //Si NO hay modo : no hacer nada
      return;
      
    int x = Math.round(ratonX);
    int y = Math.round(ratonY);
    
    if (this.modo.equals("Eliminar")) {     
      Ramal ramal = selector.seleccionarRamal(x, y, rad, circuito);
      
      Componente comp = selector.seleccionarComponente (x, y, rad, ramal);
      
      //Eliminar Componente
      if (comp != null) {
        ramal.eliminarComponente(comp);
        return;
      }
      
      //Eliminar Ramales
      if (ramal != null) {
        circuito.eliminarRamal(ramal);
        return;
      }
    }
    
    if (this.modo.equals("Seleccion")) {    //Modo Selección de Componentes
      this.ult_nodo = selector.seleccionarNodo(x, y, rad, circuito);    //Buscar nodo
      
      if (this.ult_nodo == null)                                      //Buscar Ramal
        this.ult_ramal = selector.seleccionarRamal(x, y, rad, circuito);    
      
      if (this.ult_ramal != null)
        this.ult_comp = selector.seleccionarComponente(x, y, rad, this.ult_ramal);
      
      this.modificar = !this.modificar;
      return;
    }
    
    if (this.modo.equals("Cables")) {       //Modo creación de Cables
      if (!validarColocacionCable(x, y, rad, selector, circuito)) {
        cont = 0;
        return;
      }
    
      if (this.cont == 0) {    //Primer ingreso : solo registrar valores      
        this.ult_nodo = selector.seleccionarNodo (x, y, rad, circuito);
        this.cont += 1;
        
        if (this.ult_nodo == null) {                      //Verificar si el primer click corta algún cable
          cortarCable (ratonX, ratonY, rad, circuito, selector);
          this.ult_nodo = selector.seleccionarNodo (x, y, rad, circuito);
        }
         
        //Registrar última posición presionada
        this.ult_x = x;
        this.ult_y = y;
        return;
      }
      
      Nodo nodo = selector.seleccionarNodo (x, y, rad, circuito);
      this.cont = 0;
      
      if (nodo == null) {                                 //Verificar si el último click corta algún cable
        cortarCable (ratonX, ratonY, rad, circuito, selector);
        nodo = selector.seleccionarNodo (x, y, rad, circuito);
      }
      
      if (nodo == null && this.ult_nodo == null) {        //Si no fue seleccionado ningún Nodo
        this.crearCable(x, y, this.ult_x, this.ult_y, circuito);
        return;
      }
        
      if (nodo != null && this.ult_nodo != null) {        //Si los dos extremos ya existen
        this.crearCable(nodo, this.ult_nodo, circuito);
        this.ult_nodo = null;
        return;
      }
      
      if (this.ult_nodo != null)                          //Si uno de los extremos existe
        nodo = this.ult_nodo;
        
      this.crearCable (nodo, x, y, circuito);
      return;
    }
    
    //Crear Nueva Componente : Buscar Ramal
    Ramal ramal = selector.seleccionarRamal(x, y, rad, circuito);
    
    //Validar ubicación del nuevo componente
    if (!validarColocacionComponente (x, y, rad, selector, ramal))
      return;
    
    //Obtener la Posición del Nuevo Componente
    int pos = Math.round(ramal.calcularPosComponente(x, y));
    
    //Crear el componente
    if (this.modo.equals("Resistor")) {     //Modo creación de Resistores
      this.crearResistor(pos, ramal);
      return;
    }
    
    if (this.modo.equals("Fuente")) {       //Modo creación de Fuentes
      this.crearFuente(pos, ramal);
      return;
    }
  }
  
  //---------------------------|Cambiar Orientación|---------------------------\\
  public void cambiarOrientacion () {
    this.orientacion = !(this.orientacion);
  }
  
  //---------------------------|Obtener Orientación|---------------------------\\
  public boolean obtenerOrientacion () {
    return this.orientacion;
  }
  
  //---------------------------------------------------------------------------\\
  //                               |Validaciones|                              \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Validar Corte de Cable|---------------------------\\
  public boolean validarColocacionCable (float x, float y, int rad, Selector selector, Circuito circuito) {
    //Validar si hay un cable en el camino
    Ramal ramal = selector.seleccionarRamal(x, y, rad, circuito);
    if (ramal == null)
      return true;
    
    //Validar si hay un componente en el camino
    Componente comp = selector.seleccionarComponente(x, y, rad, ramal);
    if (comp == null)    //Validar que no haya componentes en medio
      return true;
    
    return false;
  }
  
  //---------------------------|Validar Componente|---------------------------\\
  public boolean validarColocacionComponente (float x, float y, int rad, Selector selector, Ramal ramal) {
    if (ramal == null)
      return false;
    
    //Validar que no haya otro componente en el lugar
    Componente comp = selector.seleccionarComponente(x, y, rad, ramal);
    
    
    if (comp != null)
      return false;
    
     
    //Obtener la posición en el ramal asignada al Componente
    float x1 =  ramal.obtenerNodo1().obtenerPosX();
    float y1 =  ramal.obtenerNodo1().obtenerPosY();
    float x2 =  ramal.obtenerNodo2().obtenerPosX();
    float y2 =  ramal.obtenerNodo2().obtenerPosY();
      
    float[] punto = Calculadora_geo.puntoPerpendicular(x, y, x1, y1, x2, y2);
    int pos = Math.round(Calculadora_geo.calcularDistancia(punto[0], punto[1], x1, y1));
    
    //Verificar que el Componente no exceda el tamaño del Ramal
    if (ramal.longitud() <= pos + 300)
      return false;
    
    //Verificar que en la posición final del Componente no exista otro Componente
    float[] punto_final = Calculadora_geo.calcularPuntoRecta (x1, y1, x2, y2, pos + 300);
    comp = selector.seleccionarComponente(punto_final[0], punto_final[1], rad, ramal);
    
    if (comp != null)
      return false;
      
    return true;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                         |Creación de Componentes|                         \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Crear Cable|---------------------------\\
  private void crearCable (int x1, int y1, int x2, int y2, Circuito circuito) {
    if (Calculadora_geo.calcularDistancia(x1, y1, x2, y2) <= 0)
      return;
    
    //Crear Nodos
    Nodo nodo1 = new Nodo (x1, y1);
    Nodo nodo2 = new Nodo (x2, y2);
    circuito.agregarNodo(nodo1);
    circuito.agregarNodo(nodo2);
    
    //Crear Ramal
    Ramal ramal = new Ramal (nodo1, nodo2);
    circuito.agregarRamal (ramal);
  }
  
  private void crearCable (Nodo nodo, int x, int y, Circuito circuito) {
    //Crear Nodo
    Nodo nodo1 = new Nodo (x, y);
    circuito.agregarNodo(nodo1);
    
    //Crear Ramal
    Ramal ramal = new Ramal (nodo, nodo1);
    circuito.agregarRamal (ramal);
  }
  
  private void crearCable (Nodo nodo1, Nodo nodo2, Circuito circuito) {
    if (nodo1 == nodo2)  //No colocar cable de longitud 0
      return;
    
    //Crear Ramal
    Ramal ramal = new Ramal (nodo1, nodo2);
    circuito.agregarRamal (ramal);
  }
  
   //---------------------------|Separar Cable|---------------------------\\
  private void cortarCable (float x, float y, int rad, Circuito circuito, Selector selector) {
    Ramal ramal = selector.seleccionarRamal(x, y, rad, circuito);
    if (ramal == null)    //Validar que haya cable para cortar
      return;
      
    Componente comp = selector.seleccionarComponente(x, y, rad, ramal);
    if (comp != null)    //Validar que no haya componentes en medio
      return;
    
    //Obtener la Posición del Corte en el Ramal
    int pos = Math.round(ramal.calcularPosComponente(x, y));
        
    circuito.separarRamal(ramal, pos);
  }
  
  //--------------------------|Crear Resistor|-------------------------\\
  private void crearResistor (int pos, Ramal ramal) {
    ramal.agregarComponente(new Resistor(pos, this.valor));
  }
  
  //---------------------------|Crear Fuente|--------------------------\\
  private void crearFuente (int pos, Ramal ramal) {
    ramal.agregarComponente(new Fuente(pos, this.valor, this.orientacion));
  }
  
  
  //---------------------------------------------------------------------------\\
  //                          |Modificación de Modos|                          \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Activar Modo Fuente|-------------------------\\
  public void modoSeleccion (boolean activar) {
    desactivarModos();
    if (activar)
      this.modo = "Seleccion";
  }
  
  //-------------------------|Activar Modo Fuente|-------------------------\\
  public void modoCables (boolean activar) {
    desactivarModos();
    if (activar)
      this.modo = "Cables";
  }
  
  //-------------------------|Activar Modo Fuente|-------------------------\\
  public void modoResistor (boolean activar, int resistencia) {
    desactivarModos();
    if (activar) {
      this.modo = "Resistor";
      this.valor = resistencia;  
    }
  }
  
  //-------------------------|Activar Modo Fuente|-------------------------\\
  public void modoFuente (boolean activar, int voltaje) {
    desactivarModos();
    if (activar) {
      this.modo = "Fuente";
      this.valor = voltaje;
    }
  }
  
  //-------------------------|Activar Modo Fuente|-------------------------\\
  public void modoEliminacion (boolean activar) {
    desactivarModos();
    if (activar)
      this.modo = "Eliminar";
  }
  
  //-------------------------|Activar Modo Fuente|-------------------------\\
  public void desactivarModos () {
    this.modo = "Ninguno";
    this.modificar = false;
    this.valor = 0;
    this.cont = 0;
    this.ult_nodo = null;
    this.ult_ramal = null;
    this.ult_comp = null;
  }
  
  //-------------------------|Obtener Modo|-------------------------\\
  public String obtenerModo () {
    return this.modo;
  }
}
