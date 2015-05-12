#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//Value below which value is effectively 0
#define ZERO -8

//Coordinate structure
struct coord
{
    double x, y;
};

//Returns inverse of 2D matrix
vector< vector<double> > inverse(vector< vector<double> > m);

//Returns determinant of 2D matrix
double det(vector< vector<double> > m);

//Returns transpose of given 2D matrix
vector< vector<double> > transpose(vector< vector<double> > m);

//Returns matrix product of 2 matrices
vector< vector<double> > multiply(vector< vector<double> > a, vector< vector<double> > b);

//Print matrix
void print(vector< vector<double> > m);

//Row operation, interchange two rows
vector< vector<double> > interchange(vector< vector<double> > a, int row1, int row2);

//Row operation, multiply each element of row by constant
vector< vector<double> > multiplyRow(vector< vector<double> > a, int row, double constant);

//Row operation, replace row with sum of itself and a constant multiple of another row
vector< vector<double> > replaceRow(vector< vector<double> > a, int row1, int row2, double constant);

int main()
{
    int n = -1, k = -1;

	while (k != 0)
	{
		//Enter degree of polynomial
		cout << "Enter degree (0 to exit): ";
		cin >> k;

		while (k < 0)
		{
			cout << "Retry: ";
			cin >> k;
		}

		if (k == 0)
			return 0;

		//Get n
		cout << "Enter number of data points (minimum " << 2 * k << "): ";
		cin >> n;

		while (n < 2 * k)
		{
			cout << "Retry: ";
			cin >> n;
		}

		//Create Vandermonde matrix of size n x k
		vector<double> t(k + 1, 0);
		vector< vector<double> > vandermonde(n, t);

		//Create transposed matrix
		vector< vector<double> > transposed(n, t);

		//Create V^T * V
		vector< vector<double> > multiplied;

		//Create inverse
		vector< vector<double> > inversed;

		//Create Y matrix
		vector<double> te(1, 0);
		vector< vector<double> > y(n, te);

		//Create solution matrix
		vector< vector<double> > solution(k, te);

		//Get data
		vector<coord> data;
		cout << "Enter data in the form of X Y: " << endl;
		coord temp;
		for (int i = 0; i < n; i++)
		{
			cout << i+1 << ": ";
			cin >> temp.x >> temp.y;
			data.push_back(temp);
			y[i][0] = temp.y;
		}

		//Construct Vandermonde matrix
		for (int i = 0; i < n; i++)
		{
			//Increasing exponent for rest
			for (int j = 0; j <= k; j++)
			{
				vandermonde[i][j] = pow(data[i].x, k - j);
			}
		}

		//Print data
		cout << endl << "Data: " << endl;
		for (int i = 0; i < data.size(); i++)
		{
			cout << "(" << data[i].x << ", " << data[i].y << ")" << endl;
		}
		cout << endl;

		//Print Vandermonde
		cout << "V: " << endl;
		print(vandermonde);

		//Find transpose
		transposed = transpose(vandermonde);

		//Print transposed
		cout << "V^T: " << endl;
		print(transposed);

		//Find multiplied
		multiplied = multiply(transposed, vandermonde);

		//Print multiplied
		cout << "V^T * V: " << endl;
		print(multiplied);

		//Find inversed
		inversed = inverse(multiplied);

		//Print inversed
		cout << "(V^T * V)^-1: " << endl;
		print(inversed);

		//Multiply transposed by y
		y = multiply(transposed, y);

		//Print y
		cout << "V^T * Y: " << endl;
		print(y);

		//Find solution
		solution = multiply(inversed, y);

		//Print solution
		cout << "C: " << endl;
		print(solution);
	}
}

//Returns inverse of 2D matrix
vector< vector<double> > inverse(vector< vector<double> > m)
{
	//*GAUSS JORDANIAN ELIMINATION*******
	
	//Construct identity matrix
	vector< vector<double> > id = m;
	for(int i = 0; i < id.size(); i++)
	{
		for(int j = 0; j < id[i].size(); j++)
		{
			if(i==j)
				id[i][j] = 1;
			else
				id[i][j] = 0;
		}
	}

	//Solve row by row
	double multiplier = 0;
	for(int i = 0; i < id.size(); i++)
	{
		//Make location of pivot 1
		multiplier = 1 / m[i][i];
		m = multiplyRow(m, i, multiplier);
		id = multiplyRow(id, i, multiplier);

		//Make all other locations 0
		for(int j = 0; j < id.size(); j++)
		{
			if (j == i)
				continue;

			multiplier = -m[j][i];
			m = replaceRow(m, j, i, multiplier);
			id = replaceRow(id, j, i, multiplier);
		}
	}

	return id;
}

//Returns determinant of 2D matrix
double det(vector< vector<double> > m)
{
    //Base case for 2x2
    if(m.size() == 2 && m[0].size() == 2)
    {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }

	//Base case for 1x1
	else if(m.size() == 1 && m[0].size() == 1)
	{
		return m[0][0];
	}

    //Otherwise, call recursively
    else
    {
        vector< vector<double> > t;
        vector<double> temp;
        double d = 0;

        for(int i = 0; i < m[0].size(); i++)
        {
            //Delete t
            t.clear();

            //Push into t elements not in row nor column
            for(int j = 1, c = 0; j < m.size(); j++, c++)
            {
                t.push_back(temp);

                for(int k = 0; k < m[j].size(); k++)
                {
                    if(k != i)
                    {
                        t[c].push_back(m[j][k]);
                    }
                }
            }

            //Add to determinant if even
            if(i%2 == 0)
            {
                d += m[0][i] * det(t);
            }

            //Subtract if odd
            else
            {
                d -= m[0][i] * det(t);
            }
        }

        return d;
    }
}

//Returns transpose of given 2D matrix
vector< vector<double> > transpose(vector< vector<double> > m)
{
	//Create transposed with flipped sizes
	vector<double> temp(m.size());
    vector< vector<double> > t(m[0].size(), temp);

	//Assigns
    for(int i = 0; i < m.size(); i++)
    {
        for(int j = 0; j < m[i].size(); j++)
        {
            t[j][i] = m[i][j];
        }
    }

    return t;
}

//Returns matrix product of 2 matrices
vector< vector<double> > multiply(vector< vector<double> > a, vector< vector<double> > b)
{
    //Create product matrix
    vector<double> temp(b[0].size());
    vector< vector<double> > product(a.size(), temp);

	//Go through each row
    for(int i = 0; i < product.size(); i++)
    {
        for(int j = 0; j < product[i].size(); j++)
        {
            double sum = 0;

			//Adds up
            for(int k = 0; k < b.size(); k++)
            {
                sum += a[i][k] * b[k][j];
            }

            product[i][j] = sum;
        }
    }

    return product;
}

//Print matrix
void print(vector< vector<double> > m)
{
    //Go through rows
    for(int i = 0; i < m.size(); i++)
    {
        //Go through columns
        for(int j = 0; j < m[i].size(); j++)
        {
            cout << m[i][j] << " ";
        }

        cout << endl;
    }

	cout << endl;
}

//Row operation, interchange two rows
vector< vector<double> > interchange(vector< vector<double> > a, int row1, int row2)
{
	vector<double> temp = a[row1];
	a[row1] = a[row2];
	a[row2] = temp;
	return a;
}

//Row operation, multiply each element of row by constant
vector< vector<double> > multiplyRow(vector< vector<double> > a, int row, double constant)
{
	for(int i = 0; i < a[row].size(); i++)
	{
		a[row][i] *= constant;
	}

	return a;
}

//Row operation, replace row with sum of itself and a constant multiple of another row
vector< vector<double> > replaceRow(vector< vector<double> > a, int row1, int row2, double constant)
{
	for(int i = 0; i < a[row1].size(); i++)
	{
		a[row1][i] += a[row2][i] * constant;
	}

	return a;
}