/*
|====================================================================|
*                             |Selector|         
* Descripción:                                                        
*   Clase utilizada para el análisis de las colisiones en todo el
*   Circuito, arrojando los componentes con los que se encuentra
*   cerca del valor solicitado
|====================================================================|
*/

public class Selector {
  Nodo nodo_slc;            //Último Nodo seleccionado
  Ramal ramal_slc;          //Último Ramal seleccionado
  Componente comp_slc;      //Último Componente seleccionado
  
  public Selector () {
    this.nodo_slc = null;
    this.ramal_slc = null;
    this.comp_slc = null;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                   |Seleccionar Componentes del Circuito|                  \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Detectar Selección|----------------------------\\
  public void detectarSeleccion (float ratonX, float ratonY, int rad, Circuito circuito) {
    this.nodo_slc = this.seleccionarNodo(ratonX, ratonY, rad, circuito);
    this.ramal_slc = this.seleccionarRamal(ratonX, ratonY, rad, circuito);
    if (this.ramal_slc != null)
      this.comp_slc = this.seleccionarComponente(ratonX, ratonY, rad, this.ramal_slc);
  }
  
  //---------------------------|Seleccionar Nodo|----------------------------\\
  public Nodo seleccionarNodo (float px, float py, int rad, Circuito circuito) {
    ArrayList<Nodo> nodos = circuito.obtenerNodos();

    
    for (int i = 0; i < nodos.size(); i++) {
      int x = nodos.get(i).obtenerPosX();
      int y = nodos.get(i).obtenerPosY();
      if (Calculadora_geo.busqRadial(px, py, x, y, rad))
        return nodos.get(i);
    }
    
    return null;
  }
  
  //---------------------------|Seleccionar Ramal|---------------------------\\
  public Ramal seleccionarRamal (float px, float py, int rad, Circuito circuito) {
    ArrayList<Ramal> ramales = circuito.obtenerRamales();
    
    for (int i = 0; i < ramales.size(); i++) {
      int x1 = ramales.get(i).obtenerNodo1().obtenerPosX();
      int y1 = ramales.get(i).obtenerNodo1().obtenerPosY();
      int x2 = ramales.get(i).obtenerNodo2().obtenerPosX();
      int y2 = ramales.get(i).obtenerNodo2().obtenerPosY();
      if (Calculadora_geo.busqLineal(px, py, x1, y1, x2, y2, rad))
        return ramales.get(i);
    }
    
    return null;
  }
  
  //-------------------------|Seleccionar Componente|------------------------\\
  public Componente seleccionarComponente (float px, float py, int rad, Ramal ramal) {
    if (ramal == null)
      return null;
    
    ArrayList<Componente> componentes = ramal.obtenerComponentes();
    int[] x = new int[] {ramal.obtenerNodo1().obtenerPosX(), ramal.obtenerNodo2().obtenerPosX()};
    int[] y = new int[] {ramal.obtenerNodo1().obtenerPosY(), ramal.obtenerNodo2().obtenerPosY()};
    
    for (int i = 0; i < componentes.size(); i++) {
      int[][] puntos = Calculadora_geo.aproxMatriz(componentes.get(i).obtenerPuntos(x, y));
      if (Calculadora_geo.busqLineal(px, py, puntos[0][0], puntos[1][0], puntos[0][1], puntos[1][1], rad))
        return componentes.get(i);
    }
    
    return null;
  }
}
