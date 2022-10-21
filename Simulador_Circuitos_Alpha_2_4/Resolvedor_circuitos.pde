/*
|====================================================================|
*                      |Resolvedor de Circuitos|     
*                          (Clase estática)
* Descripción:                                                        
*   Clase utilizada para 
|====================================================================|
*/

public static class Resolvedor_circuitos {
  
  private static int num_mallas;
  
  
  //-------------------------|Resolver Circuito|-------------------------\\
  public static void resolverCircuito (Circuito circuito) {
    ArrayList<Circuito> mallas = obtenerMallas(circuito);  //Obtener las Corrientes
    float[] corrientes = obtenerCorrientes(mallas);
    
    //Asignar las corrientes
    for (int i = 0; i < mallas.size(); i++) {
      ArrayList<Ramal> ramales = mallas.get(i).obtenerRamales();
      
      for (Ramal ramal : ramales)
        ramal.modificarCorriente(corrientes[i]);    //Sumar corriente a cada ramal
    }
    
    //Ajustar las corrientes a positivo y su dirección
    for (int i = 0; i < mallas.size(); i++) {
      ArrayList<Ramal> ramales = mallas.get(i).obtenerRamales();
      
      for (Ramal ramal : ramales)
        ramal.ajustarDirCorriente();    //Sumar corriente a cada ramal
    }
  }
  
  //-------------------------|Obtener Mallas del Circuito|-------------------------\\
  public static ArrayList<Circuito> obtenerMallas (Circuito circuito) {
    Vector<Integer>[] matriz_ady = obtenerMatrizAdy (circuito);
    Vector<Integer>[] ciclos = new Vector[50];
    num_mallas = 0;
    
    int[] col = new int[50];
    int[] par = new int[50];
    
    dfs_cycle(matriz_ady, ciclos, 1, 0, col, par);
    
    println(num_mallas);
    for (Vector<Integer> ciclo : ciclos) {
      if (ciclo == null)
        continue;
      
      println("Ciclo: " + ciclo);
    }
    
    return null;
  }
  
  
  
  //-------------------------|Generar Matriz de Adyacencia|-------------------------\\
  public static Vector<Integer>[] obtenerMatrizAdy (Circuito circuito) {
    ArrayList<Nodo> nodos = circuito.obtenerNodos();
    Vector<Integer>[] matriz = new Vector[nodos.size()];
    
    for (int i = 0; i < nodos.size(); i++) {
      ArrayList<Nodo> vecinos = nodos.get(i).obtenerVecinos();
      matriz[i] = new Vector<Integer>();
      
      for (Nodo vecino : vecinos)  //Agregar vecinos a matriz de adyacencia
        matriz[i].add(nodos.indexOf(vecino));
    }
    
    for (Vector<Integer> i : matriz)
      println(i);
    
    return matriz;
  }
  
  //-------------------------|Obtener Ciclos del Circuito|-------------------------\\
  public static void dfs_cycle(Vector<Integer>[] graph, Vector<Integer>[] cycles, int u, int p, int[] col,int[] par) {
    // already (completely) visited vertex.
    if (col[u] == 2)
      return;
 
        // seen vertex, but was not completely visited -> cycle detected.
        // backtrack based on parents to find the complete cycle.
        if (col[u] == 1)
        {
 
             
              Vector<Integer> v = new Vector<Integer>();
            int cur = p;
              v.add(cur);
 
            // backtrack the vertex which are
            // in the current cycle thats found
            while (cur != u)
            {
                cur = par[cur];
                v.add(cur);
            }
              cycles[num_mallas] = v;
            num_mallas++;
            return;
        }
        par[u] = p;
 
        // partially visited.
        col[u] = 1;
 
        // simple dfs on graph
        for (int v : graph[u])
        {
 
            // if it has not been visited previously
            if (v == par[u])
            {
                continue;
            }
            dfs_cycle(graph, cycles, v, u, col, par);
        }
 
        // completely visited.
        col[u] = 2;
  }
  
  //---------------------------------------------------------------------------\\
  //                     |Resolución Sistema de Ecuaciones|                    \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Obtener Corrientes de Malla|-------------------------\\
  public static float[] obtenerCorrientes (ArrayList<Circuito> mallas) {
      double[][] matriz = generarMatrizMalla(mallas);
      double[][] coeficientes = new double[mallas.size()][mallas.size()];
      
      for (int i = 0; i < matriz.length - 1; i++)
        coeficientes[i] = matriz[i];
      
      double[] valores = matriz[matriz.length - 1];
        
      double[] res = resolverSistema(coeficientes, valores);
      float[] corrientes = new float[res.length];
        
      for(int i = 0; i < res.length; i++) {
        corrientes[i] = (float)res[i];
        System.out.print(res[i]+" ");
      }
      
      return corrientes;
  }
  
  //-------------------------|Resolver Sistema de Ecuaciones|-------------------------\\
  public static double[] resolverSistema(double[][] matrix, double[] coef) {
    try {
      RealMatrix coefficients = new Array2DRowRealMatrix(matrix ,false);
      DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
      RealVector constants = new ArrayRealVector(coef, false);
      RealVector solution = solver.solve(constants);  
      return solution.toArray();
    } catch (SingularMatrixException e) {
      return new double[0];
    }
  }
  
  
  //---------------------------------------------------------------------------\\
  //                     |Generación Matriz de Ecuaciones|                     \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Generar Matriz de Ecuaciones|-------------------------\\
  public static double[][] generarMatrizMalla (ArrayList<Circuito> mallas) {
    double[][] ecuaciones = new double[mallas.size() + 1][];
    ecuaciones[mallas.size()] = new double [mallas.size()];
    
    for (int i = 0; i < mallas.size(); i++) {   //Generar la ecuación de cada malla
      double[] ecua = generarEcuacionMalla(mallas, i);
      double[] ecua_cor = new double [ecua.length - 1];
      
      for (int j = 0; j < ecua.length - 1; j++) {
        ecua_cor[j] = ecua[j];
      }
      
      ecuaciones[i] = ecua_cor;    //Guardar coeficientes
      ecuaciones[mallas.size()][i] = ecua[ecua.length - 1];  //Guardar valores - última fila
    }
      
    return ecuaciones;
  }
  
  //-------------------------|Generar Ecuación de Malla|-------------------------\\
  public static double[]generarEcuacionMalla (ArrayList<Circuito> mallas, int pos) {
    double[] ecuacion = new double [mallas.size() + 1];
    
    ArrayList<Ramal> ramales = mallas.get(pos).obtenerRamales();
    
    for (Ramal ramal : ramales) {
      //Agregar Componentes del Ramal a la Ecuación
      ecuacion[pos] = ecuacion[pos] + ramal.obtenerResistencia();
      ecuacion[mallas.size()] = ecuacion[mallas.size()] + ramal.obtenerFEM();
      
      for (int i = 0; i < mallas.size(); i++) {    //Verificar si el ramal hace parte de otra malla
         if (mallas.get(i) == mallas.get(pos))    //No realizar operación en el mismo circuito
           continue;
           
         if (ramalEnCircuito(mallas.get(i), ramal))  //Si lo es, incluirlo en la ecuación
           ecuacion[i] = ecuacion[i] - ramal.obtenerResistencia();
      }
    }
    
    return ecuacion;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                               |Validaciones|                              \\
  //---------------------------------------------------------------------------\\
  
  //-------------------------|Verificar si un Componente está en el Circuito|-------------------------\\
  public static boolean componenteEnCircuito (Circuito circuito, Componente componente) {
    ArrayList<Ramal> ramales = circuito.obtenerRamales();
    
    for (Ramal ramal : ramales) {      //Verificar Cada Ramal del Circuito
      ArrayList<Componente> componentes = ramal.obtenerComponentes();
      
      if (componentes.contains(componente))  //Verificar los componentes de cada ramal
        return true;
    }
      
    return false;
  }
  
  //-------------------------|Verificar si un Ramal está en el Circuito|-------------------------\\
  public static boolean ramalEnCircuito (Circuito circuito, Ramal ramal) {
    return circuito.obtenerRamales().contains(ramal);
  }
}
