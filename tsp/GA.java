//	遺伝的アルゴリズム(巡回セールスマン問題)　　本体
package tsp;

import java.awt.Point;
import java.util.*;

//	遺伝的アルゴリズム　（巡回セールスマン問題）
public class GA {
	enum		CType  { Cross, Cut };
	enum		MType  { Change, Twist };
	private	CType	ctype = CType.Cut;
	private	MType	mtype = MType.Twist;
	
	private	int					generation;	// 世代数
	private	ArrayList<Point>	original;		// 元の座標列
	private	ArrayList<Point>	bestpath;		// 最善パス
	private	int					scale;			// ±scale
	private	double				epsilon = scale * 0.000001;
	private	int					pathlength;		// 点の数
	private	ArrayList<Individual>		pool;	// 個体プール
	private	int					poolinitsize = 50;
	private	int					selectsize = 10;
	private	int					crosssize  = 30;
	private	int					mutatesize = 10;
	
	private	LineEq				cutline = null;
	
	// コンストラクタ：点列を与えて初期化
	GA( ArrayList<Point> ps,int sc ) {
		original = ps;		// 位置データのセット
		scale = sc;
		pathlength = original.size();
		generation = -1;		// GA生成時は、世代を -1 にしておく
		bestpath = null;
	}

	// getメソッド　　　点列、道程、世代数　を得る
	ArrayList<Point>	getPath()	{	return	bestpath;		}
	int			getMinLength() 	{	return	(int)(pool.get(0).distance);	}
	int			getGeneration() 	{	return	generation;				}
		
	// 進化を進める
	void	stepforward() {
		if( generation < 0 )		// プールの初期化（第0世代）
			initializepool();				
		onegeneration();			// 進化を1世代進める	
		++generation;
		bestpath = pool.get(0).journey;
	}
	
	// 1世代分の処理
	private void	onegeneration() {
		ArrayList<Individual>	nextpool = new ArrayList<Individual>();
		
		for( int i=0; i<selectsize; ) {	// 選択
			Individual indv = pool.get(i);
			if( find(indv,nextpool) >= 0  )
				continue;				
			nextpool.add(indv); ++i;
		}

		for( int i=0; i<crosssize; ) {		// 交叉
			Individual indv0 = selectRandomfromPool();
			Individual indv1 = selectRandomfromPool();
			if( indv0 == indv1 )
				continue;
			Individual indv = genCross(indv0,indv1);
			if( find(indv,nextpool) >= 0 )
				continue;
			nextpool.add(indv); ++i;
		}

		for( int i=0; i<mutatesize; ) {		// 突然変異
			Individual indv = genMutation(selectRandomfromPool());
			if( find(indv,nextpool) >= 0 )
				continue;	
			nextpool.add(indv); ++i;
		}

		pool = nextpool;
		sortpool();
	}
	
	// 染色体が個体プールに含まれるかの検査   indexを返す
	private int find( Individual indv, ArrayList<Individual> list ) {
		int	idx = 0;
		for( Individual id : list ) {
			if( id.equals(indv) )
				return	idx;
			++idx;
		}
		return	-1;
	}
	
	//	プールの初期化
	void	initializepool() {
		pool = new ArrayList<Individual>();
		for( int i=0; i<poolinitsize; ++i ) {
			Individual id = new Individual();
			pool.add(id);
		}
		sortpool();
		generation = 0;
		bestpath = pool.get(0).journey;
	}

	//  交叉染色体の生成
	private Individual genCross( Individual indv0, Individual indv1 ) {
		switch( ctype ) {
		case Cross:
			return	genCrossOrder(indv0,indv1);
		case Cut:
			return	genCrossCutting(indv0,indv1);
		}
		return	null;		// dummy
	}

	//	２点で交叉する染色体の生成　　g0-g1-g0	【タイプ１　２点交差】
	private	Individual genCrossOrder( Individual indv0, Individual indv1 ) {
		int	sp[], xlen;
		do {
			sp = randomposition2(pathlength);
			xlen = sp[1]-sp[0];
		} while( xlen < pathlength/4 || xlen > pathlength*3/4 );
		int[]	oc = Arrays.copyOf(indv0.ordercode,indv0.ordercode.length);

		for( int i=sp[0]; i<sp[1]; ++i )
			oc[i] = indv1.ordercode[i];
		
		return		new Individual(oc);
	}
	
	// 直線で盤面分割する方法				【タイプ２　直線分割】
	private	Individual genCrossCutting( Individual indv0, Individual indv1 ) {
		// 切断する直線を決める
		do {
			double	th = 2*Math.PI*Math.random();
			double  a = Math.cos(th);
			double  b = Math.sin(th);		
			double	x = scale*(Math.random() - 0.5);
			double	y = scale*(Math.random() - 0.5);
			double  c = -(a*x+b*y);
			cutline = new LineEq(a,b,c);
		} while( !dividecheck(cutline) );
		
		// ２つの染色体を直線で切り、残す側のセグメントをまとめる
		ArrayList<ArrayList<Point>>	segs = segmentalize(indv0,cutline,1);
		ArrayList<ArrayList<Point>>	segs1 = segmentalize(indv1,cutline,-1);	
		for( ArrayList<Point> seg : segs1 )	// Join
			segs.add(seg);
		
		// セグメントの端点を集める
		ArrayList<Point>	endpoints = new 	ArrayList<Point>();
		for( ArrayList<Point> seg : segs ) {
			Point sp = seg.get(0);
			Point ep = seg.get(seg.size()-1);
			endpoints.add(sp);
			if( !sp.equals(ep) )
				endpoints.add(ep);
		}
		
		// セグメントを繋ぐ
		ArrayList<Point>	newpath = new ArrayList<Point>();
		Point	np = endpoints.get(0);		// 端点
		while( segs.size() > 0) {			// セグメントが尽きるまで繰り返す
			// 	端点np を含むセグメントを探す
			ArrayList<Point> nextseg = null;
			for( ArrayList<Point> seg : segs ) {
				if( seg.get(0).equals(np) || seg.get(seg.size()-1).equals(np) ) {
					nextseg = seg;
					segs.remove(nextseg);
					break;
				}
			}
			// 見つけたセグメントを、newpathに繋ぐ
			if( nextseg.get(0).equals(np) ) {
				for( Point p : nextseg )
					newpath.add(p);
			} else {
				while( nextseg.size()>0 ) {
					newpath.add(nextseg.get(nextseg.size()-1));
					nextseg.remove(nextseg.size()-1);
				}
			}
			// 端点のリストから、不要になった端点を消す
			endpoints.remove(np);
			if( ! newpath.get(newpath.size()-1).equals(np)) {
				np = newpath.get(newpath.size()-1);
				endpoints.remove(np);
			}
			// 新たに繋ぐ端点npを求める
			Point	pathend = newpath.get(newpath.size()-1);
			double	mindistance = Double.MAX_VALUE;		// 十分大きな値
			for( Point p : endpoints ) {
				double dist = p.distance(pathend);
				if( dist < mindistance ) {
					dist = mindistance;
					np = p;
				}
			}
		}
		return	new Individual(newpath);
	}
	
	private boolean	dividecheck( LineEq line ) {
		int		mincount = pathlength / 10 + 1;
		int		pside=0, mside=0;
		for( Point p : original ) {
			if( line.whichside(p)>0 )
				++pside;
			else
				++mside;
			if( pside>mincount && mside>mincount )
				return	true;
		}
		return	false;
	}

	// 経路を切断線で切り、side側の経路断片のリストを返す。
	private ArrayList<ArrayList<Point>> segmentalize( Individual indv, LineEq line, int side ) {
		ArrayList<ArrayList<Point>>	segments = new ArrayList<ArrayList<Point>>();
		ArrayList<Point>	seg = null;
		
		for( Point p : indv.journey ) {
			if( line.whichside(p) != side ) {
				if( seg != null )
					segments.add(seg);
				seg = null;
				continue;
			}
			if( seg == null )
				seg = new ArrayList<Point>();
			seg.add(p);
		}
		if( seg != null )
			segments.add(seg);
		
		return	segments;
	}
	
	//	突然変異
	private Individual	genMutation( Individual indv ) {
		switch( mtype ) {
		case Change:
			return	getMutationChange(indv);
		case Twist:
			return	getMutationTwist(indv);
		}
		return	null;
	} 
	// 順序表現の2箇所を変更する
	private Individual	getMutationChange( Individual indv ) {
		int[]	oc = Arrays.copyOf(indv.ordercode,indv.ordercode.length);
	
		for( int i=0; i<2; ++i ) {
			int	xp = (int) (Math.random()*(pathlength));
			int	m  = (int) (Math.random()*(pathlength-xp));	
			oc[xp] = m;
		}
		
		return		new Individual(oc);
	}
	
	//	ランダムな１区間内を逆順に並べる
	private Individual	getMutationTwist( Individual indv ) {
		ArrayList<Point>	path = new ArrayList<Point>();
		int	sp[] = randomposition2(pathlength);
	
		for( int i=0; i<sp[0]; ++i )
			path.add(indv.journey.get(i));
		for( int i=sp[0]; i<sp[1]; ++i )
			path.add(indv.journey.get(sp[0]+sp[1]-1-i));
		for( int i=sp[1]; i<pathlength; ++i )
			path.add(indv.journey.get(i));

		return		new Individual(path);
	}
	
	// プール内を、パス長の順にソート
	private void	sortpool() {
		Collections.sort(pool,new IndividualComparator());
	}
	
	// ２つのランダムな異なる添字を戻す
	private int[]	randomposition2( int len ) {
		int[]	rx = new int[2];
		do {
			rx[0] = (int) (Math.random()*(len));
			rx[1] = (int) (Math.random()*(len));
		} while( rx[0] >= rx[1] );
		return	rx;
	}

	// 個体プールからランダムに選択する　（優先度あり）
	private Individual	selectRandomfromPool() {
		int	idx = (int)(Math.random()*pool.size()/2);
		return	pool.get(idx);
	}
	
	// 個体（染色体） =============================================================
	class  Individual {						// 染色体
		ArrayList<Point>	journey;		// 　パス表現(点列)
		int[]				ordercode;		// 　順序表現
		double				distance;		// 道程

		//	コンストラクタ
		Individual() 							{	makeGene( makeRandomPath() );	}
		Individual( ArrayList<Point> list ) 	{ 	makeGene( list );		}	
		Individual( int[] oc ) {					// 順序表現から個体を生成
			ordercode = oc;
			order2journey();
			setDistance();
		}
		private void makeGene( ArrayList<Point> list ) {	// 点列から個体を生成
			journey = list;
			if( ctype == CType.Cross || mtype == MType.Change )
				journey2order();
			setDistance();			
		}
		
		// ランダムな経路を作る
		private	ArrayList<Point>	makeRandomPath() {
			ArrayList<Point>	work = (ArrayList<Point>) original.clone();	
			ArrayList<Point>  points = new ArrayList<Point>();
			while( work.size() > 0 ) {
				int	idx = (int) (Math.random()*work.size());
				points.add(work.get(idx));
				work.remove(idx);
			}
			return	points;
		}
		
		// 経路の全長を求める
		void	setDistance() {
			double len = 0.0;
			if( journey == null || journey.size()<2 )
				return;
			Point	prev = journey.get(journey.size()-1);
			for( Point next : journey )	{
				len += prev.distance(next);
				prev = next;
			}
			distance = len;
		}
		
		//	順序表現からパス表現に直す
		void	order2journey() {
			ArrayList<Point>	work = (ArrayList<Point>) original.clone();
			journey = new ArrayList<Point>();
			for( int i : ordercode ) {
				journey.add(work.get(i));
				work.remove(i);
			}
		}
		
		//	パス表現から順序表現に直す
		void	journey2order() {
			ArrayList<Point>	work = (ArrayList<Point>) original.clone();
			ordercode = new int[journey.size()];
			for( int i=0; i<ordercode.length; ++i ) {
				Point	p = journey.get(i);
				ordercode[i] = work.indexOf(p);
				work.remove(p);
			}
		}
		
		//	個体の一致検査
		boolean	equals(Individual indv) {
			if( Math.abs(indv.distance - distance) > epsilon )
				return	false;
			/*		道程が非常に近いだけで判断して十分のようである */
			for( int i=0; i<pathlength; ++i )
				if( journey.get(i) != indv.journey.get(i) ) {
					return	false;
				}
			/**/
			return	true;
		}
	}
	
	// ２個体のコンパレータ
	class IndividualComparator implements Comparator {
		public int compare(Object p, Object q) {
			double	d = ((Individual)p).distance - ((Individual)q).distance;
			return d > 0.0 ? 1 : d < 0.0 ? -1 : 0;
		}
	}
}

// 切断線 ===================================================================
class	LineEq {
	double	a, b, c;
	
	LineEq(double aa, double bb, double cc ) {
		a = aa;
		b = bb;
		c = cc;
	}
	
	int		whichside( Point p ) {
		double	d = a*p.x + b*p.y + c;
		return	d >0 ? 1 : -1;
	}
}