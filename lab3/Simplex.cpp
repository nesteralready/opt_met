//
// Created by Ilya Nesterenko on 19.04.2021.
//
#include <iostream>
#include <vector>

class Simplex {
public:
  std::vector<std::vector<double>> table;
  int m, n;

  std::vector<int> basis;

  Simplex(std::vector<std::vector<double>> &source) {
    m = source.size();
    n = source[0].size();
    table.resize(m, std::vector<double>(n + m - 1,0));
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < table.size(); j++)
      {
        if (j < n)
          table[i][j] = source[i][j];
        else
          table[i][j] = 0;
      }
      if ((n + i) < table[0].size())
      {
        table[i][n + i] = 1;
        basis.push_back(n + i);
      }
    }

    n = table[0].size();

  }

  [[nodiscard]] std::vector<std::vector<double>> Calculate(std::vector<double> &result) {
    int mainCol, mainRow; //ведущие столбец и строка

    while (!IsItEnd()) {
      mainCol = findMainCol();
      mainRow = findMainRow(mainCol);
      basis[mainRow] = mainCol;

      std::vector<std::vector<double>> new_table (m, std::vector<double>(n));

      for (int j = 0; j < n; j++) {
        new_table[mainRow][j] = table[mainRow][j] / table[mainRow][mainCol];
      }

      for (int i = 0; i < m; i++) {
        if (i == mainRow) {
          continue;
        }

        for (int j = 0; j < n; j++) {
          new_table[i][j] = table[i][j] - table[i][mainCol] * new_table[mainRow][j];
        }
      }
      table = new_table;
    }

    for (int i = 0; i < result.size(); i++)
    {
      int k = getElem(i+1,basis);
      if (k != -1)
        result[i] = table[k][0];
      else
        result[i] = 0;
    }

    return table;
  }

  int getElem(int value, std::vector<int> source){
      for( int i =0; i< source.size(); i++){
          if(value == source[i]) {
            return i;
          }
      }
      return -1;
  }

bool IsItEnd() {
    bool flag = true;

    for (int j = 1; j < n; j++) {
      if (table[m - 1][j] < 0){
        flag = false;
        break;
      }
    }

    return flag;
  }

int findMainCol()
  {
    int mainCol = 1;

    for (int j = 2; j < n; j++)
      if (table[m - 1][j] < table[m - 1][mainCol])
        mainCol = j;

    return mainCol;
  }

int findMainRow(int mainCol)
  {
    int mainRow = 0;

    for (int i = 0; i < m - 1; i++)
      if (table[i][mainCol] > 0)
      {
        mainRow = i;
        break;
      }

    for (int i = mainRow + 1; i < m - 1; i++) {
      if ((table[i][mainCol] > 0)
          && ((table[i][0] / table[i][mainCol]) < (table[mainRow][0] / table[mainRow][mainCol]))) {
           mainRow = i;
      }
    }

    return mainRow;
  }


};