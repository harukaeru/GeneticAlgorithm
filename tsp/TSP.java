//	遺伝的アルゴリズム(巡回セールスマン問題)　　GUI部
package tsp;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Point2D;
import java.util.*;
import javax.swing.*;

public class TSP extends JFrame{
	Board board = new Board();
	int			count 	= 100;			// 総点数
	JTextField	countfield = new JTextField(""+count,4);
	int		loopcount = 1;			// ステップ用ループ数	
	JLabel	genlabel = new JLabel("generation");
	JLabel	lenlabel = new JLabel("length");

	public static void main(String[] args) {
		TSP frame = new TSP();
		frame.repaint();
	}

	public TSP() {
		Container cont = this.getContentPane();	
		cont.setLayout(new BorderLayout());

		cont.add(makeControlPanel(),"North");	// 制御パネル
		cont.add(board,"Center");

		setTitle("遺伝的アルゴリズムで解く巡回セールスマン問題");
		setSize(500 ,550);
		setVisible(true);	
	}

	// 制御パネル
	private	JPanel	makeControlPanel() {
		JPanel	panel = new JPanel();
		
		panel.add(countfield);
		
		JButton	pointbtn = new JButton("点");
		pointbtn.addActionListener(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						board.makenewpoints();
					}});
		panel.add(pointbtn);
		
		JButton	initbtn = new JButton("初");
		initbtn.addActionListener(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						board.init();
					}});
		panel.add(initbtn);
			
		String[]	loopcountnames = {"1", "10", "100", "1000","10000", "100000" };
		JComboBox	loopcountcombo = new JComboBox(loopcountnames);
		loopcountcombo.addActionListener(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						JComboBox cb = (JComboBox)e.getSource();
						loopcount = Integer.parseInt((String)(cb.getSelectedItem()));
					}});
		panel.add(loopcountcombo);
		
		JButton	stepbtn = new JButton("ステップ");
		stepbtn.addActionListener(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						board.stepforward(loopcount);
					}});
		panel.add(stepbtn);
		
		panel.add(genlabel);
		panel.add(new JLabel("/"));
		panel.add(lenlabel);
		
		return	panel;
	}

	//	巡回路の表示ボード
	class Board extends Canvas {
		Point			center;			// キャンバス中央の座標
		final	int		scale 	= 200;		// 領域 正方形　±200
		int				gap;				// 最短間隔	
		ArrayList<Point> positions;		// 点列
		GA				ga;					// 遺伝的アルゴリズム・オブジェクト
		
		// 描画
		public void paint(Graphics g) {
			Dimension	sz	= getSize();
			center	= new Point( sz.width/2, sz.height/2 );
			g.setColor(Color.white);
			g.fillRect(0,0,sz.width,sz.height);
			
			draw(g);
			if( ga!=null )
				genlabel.setText(""+ga.getGeneration());
		}

		private void draw(Graphics g) {
			g.setColor(Color.red);
			drawPoints(g);
			if( ga!=null )
				drawLines(g,ga.getPath());
		}
		
		private void	drawPoints(Graphics g) {
			if( positions == null )
				return;
			for( Point p : positions ) {
				int	u = center.x + p.x;
				int	v = center.y + p.y;
				g.fillRect(u-2,v-2,5,5);
			}
		}

		private void	drawLines( Graphics g, ArrayList<Point> ps) {
			if( ps == null )
				return;
			Point	p = ps.get(ps.size()-1);
			for( Point q : ps ) {
				int	px = center.x + p.x;
				int	py = center.y + p.y;	
				int	qx = center.x + q.x;
				int	qy = center.y + q.y;	
				g.drawLine(px,py,qx,qy);			
				p=q;
			}
		}	
		
		// 点列を生成し、GAオブジェクトを生成する
		public void	makenewpoints() {
			count = Integer.parseInt(countfield.getText());	// 点数
			gap = (int)Math.sqrt(scale*scale*2/count);	// 点数に依存した最短距離

			positions = new ArrayList<Point>();
			OUTLOOP: while( positions.size() < count ) {
				int	x = (int)(scale*2*Math.random()) - scale;
				int	y = (int)(scale*2*Math.random()) - scale;			
				Point	p = new Point(x,y);
				for( Point q : positions ) {
					if( p.distance(q) < gap )
						continue OUTLOOP;
				}
				positions.add(p);
			}

			ga = new GA(positions,scale);
			repaint();
		}

		// 同一点で再度実行可能にする。第0世代(初期プール設定)
		public	void	init() {
			if( ga == null )
				return;
			ga.initializepool();
			repaint();
			lenlabel.setText(""+ga.getMinLength());
		}

		// 進化をnステップ進める
		public	void stepforward(int n) {
			for( int i=0; i<n; ++i)
				ga.stepforward();

			repaint();	
			lenlabel.setText(""+ga.getMinLength());		
		}
	}
}