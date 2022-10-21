/*
|====================================================================|
*                      |Calculadora Geométrica|     
*                          (Clase estática)
* Descripción:                                                        
*   Clase utilizada para los cálculos geométricos y análisis de
*   colisiones necesarios para el funcionamiento de la interfaz
|====================================================================|
*/

public static class Calculadora_geo {
  //---------------------------------------------------------------------------\\
  //                           |Análisis Geométrico|                           \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Búsqueda Radial|----------------------------\\
  public  static boolean busqRadial (float px, float py, int x, int y, int rad) {
    float dist = sqrt((px - x)*(px - x) + (py - y)*(py - y));
    
    if (dist < rad)
      return true;
    
    return false;
  }
  
  //---------------------------|Búsqueda Lineal|----------------------------\\
  public static boolean busqLineal (float px, float py, int x, int y, int x2, int y2, int rad) {
    int[][] area_cuadrada = calcularAreaCuadrada (x, y, x2, y2, rad);
    
    return busqRect(px, py, area_cuadrada[0], area_cuadrada[1]);
  }
  
  //-------------------------|Búsqueda Rectangular|-------------------------\\
  public static boolean busqRect (float px, float py, int[] x, int[] y) {
    //Generar el área de cada triángulo formado entre cada lado y el punto
    float at1 = valorAbs(areaTriangular (px, py, x[0], y[0], x[1], y[1]));
    float at2 = valorAbs(areaTriangular (px, py, x[1], y[1], x[2], y[2]));
    float at3 = areaTriangular (px, py, x[2], y[2], x[3], y[3]);
    float at4 = areaTriangular (px, py, x[3], y[3], x[0], y[0]);
    
    //Calcular el area del rectángulo
    float area = areaTriangular(x[0], y[0], x[1], y[1], x[2], y[2]) + areaTriangular(x[0], y[0], x[2], y[2], x[3], y[3]);
    
    //Si el área de los triángulos es mayor está fuera del rectángulo
    if (Math.round(at1 + at2 + at3 + at4) <= area)  
      return true;
    
    return false;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                           |Calculo Geométrico|                            \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Calcular Area de Triángulo|----------------------------\\
  public static float areaTriangular (float x1, float y1, float x2, float y2, float x3, float y3) {
    return valorAbs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
  }
  
  //---------------------------|Calcular Pendiente de una Línea|----------------------------\\
  public static float calcularPendiente (float x1, float y1, float x2, float y2) {
    float m = (y2 - y1) / (x2 - x1);
    m = m * (-1.0);  //Correción (Processing trabaja en el 4to cuadrante)
    return m;
  }
  
  //---------------------------|Calcular Pendiente de una Línea|----------------------------\\
  public static float calcularDistancia (float x1, float y1, float x2, float y2) {
    return sqrt((y2 - y1)*(y2 - y1) + (x2 - x1)*(x2 - x1));
  }
  
  public static float[] puntoPerpendicular (float px, float py, float x1, float y1, float x2, float y2) {
    if (x1 == x2)    //Si la recta es paralela al eje x
      return new float[] {x1, py};
    
    float m = calcularPendiente(x1, y1, x2, y2);
    
    if (m == 0)      //Si la recta es paralela al eje y
      return new float[] {px, y1};
    
    //Calcular la distancia entre el punto y la recta
    float A = m;
    float B = -1;
    float C = (-y2) - (m * x2);
    
    float hipo = Calculadora_geo.calcularDistancia(px, py, x1, y1);                //Conseguir hipotenusa (distancia del ratón al nodo1)
    float cat = (Calculadora_geo.valorAbs(A*px + B*-py + C)) / (sqrt(A*A + B*B));  //Conseguir distancia mínima del ratón a la recta
    float pos = sqrt(hipo*hipo - cat*cat);                                         //Conseguir el cateto (distancia del nodo al punto más cercano al ratón)
    
    return Calculadora_geo.calcularPuntoRecta(x1, y1, x2, y2, pos);
  }
  
  //---------------------------|Calcular Punto en una Recta dada una Longitud|----------------------------\\
  public static float[] calcularPuntoRecta (float x1, float y1, float x2, float y2, float dist) {
    float longitud = calcularDistancia(x1, y1, x2, y2);
    
    if (longitud < dist)  //Si el punto está por fuera del segmento
      return new float[] {0.0, 0.0};
      
    float tx = (dist * (x2 - x1) / longitud);
    float ty = (dist * (y2 - y1) / longitud);
    
    return new float[] {x1 + tx, y1 + ty};
  }
  
  //---------------------------|Calcular puntos de un cuadrado|----------------------------\\
  //Dada una determinada línea esta función arrojará los puntos que equivalen a un Rectángulo formado
  //Al extender la línea sobre otro eje en un valor de rad
  //Esta función devuelve una matriz:
  //resultado[0] son los valores de x en los 4 puntos y resultado[1] son los valores de y
  public static int[][] calcularAreaCuadrada (int x1, int y1, int x2, int y2, int rad) {
    int[][] resultado = new int[2][];
    
    //Si la linea es paralela al eje y
    if (x1 == x2) {
      resultado[0] = new int[] {x1 - rad, x1 + rad, x2 + rad, x2 - rad};
      resultado[1] = new int[] {y1, y1, y2, y2};
      return resultado;
    }
    
    float m = calcularPendiente (x1, y1, x2, y2);
    
    //Si la línea es paralela al eje x
    if (m == 0) {
      resultado[0] = new int[] {x1, x1, x2, x2};
      resultado[1] = new int[] {y1 - rad, y1 + rad, y2 + rad, y2 - rad};
      return resultado;
    }
    
    //Calcular el área cuadrada según la pendiente
    float valor = (sqrt(1 + m*m))/(1 + m*m);
    int ty = Math.round(rad * valor);
    int tx = Math.round(rad * valorAbs(m) * valor);
    
    //Si la línea está inclinada
    resultado[0] = new int[] {x1 - tx, x1 + tx, x2 + tx, x2 - tx};
    
    if (m > 0) {    //Inclinada hacia la derecha
      resultado[1] = new int[] {y1 -  ty, y1 + ty, y2 + ty, y2 - ty};
    } else {        //Inclinad hacia la izquierda
      resultado[1] = new int[] {y1 + ty, y1 - ty, y2 - ty, y2 + ty};
    }
    
    return resultado;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                           |Correción de Datos|                            \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Valor Absoluto|----------------------------\\
  public static int valorAbs (int num) {
    if (num >= 0)
      return num;
    return num * -1;
  }  
  
  //---------------------------|Valor Absoluto|----------------------------\\
  public static float valorAbs (float num) {
    if (num >= 0)
      return num;
    return num * -1;
  }  
  
  //---------------------------|Aproximar Matriz|----------------------------\\
  public static int[][] aproxMatriz (float[][] mat) {
    int[][] nuevo = new int [mat.length][mat[0].length];
    
    for (int i = 0; i < mat.length; i++) {
      for (int j = 0; j < mat[i].length; j++) {
        nuevo[i][j] = Math.round(mat[i][j]);
      }
    }
    
    return nuevo;
  }
  
  //---------------------------|Aproximar Matriz|----------------------------\\
  public static int[] aproxVector (float[] vec) {
    int[] nuevo = new int [vec.length];
    
    for (int i = 0; i < vec.length; i++)
      nuevo[i] = Math.round(vec[i]);
    
    return nuevo;
  }
}
