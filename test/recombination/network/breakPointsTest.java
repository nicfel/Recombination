package recombination.network;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class breakPointsTest {
	

    @Test
    public void combinedTest() {
    	BreakPoints breakPoints = new BreakPoints(1000);
    	
    	Assert.assertEquals(breakPoints.toString(), "0-999");
    	
    	breakPoints.computeLeftAndRight(500);   
    	
    	BreakPoints left = breakPoints.getLeft();
    	BreakPoints right = breakPoints.getRight();
    	
    	
    	
    	Assert.assertEquals(left.toString(), "0-500");
    	Assert.assertEquals(right.toString(), "501-999");

    	BreakPoints bp = breakPoints.getLeft().copy();
    	bp.computeLeftAndRight(250);
    	bp = bp.getLeft();
    	bp.computeLeftAndRight(100);
    	
    	
    	bp = bp.getRight();

    	List<Integer> bp2_list = new ArrayList<>();
    	bp2_list.add(150);
    	bp2_list.add(350);
    	
    	bp2_list.add(450);
    	bp2_list.add(550);
    	    	
    	BreakPoints bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	
    	bp.or(bp2);
    	
    	Assert.assertEquals(bp.toString(), "101-350,450-550");
    	
    	bp.computeLeftAndRight(400);
    	Assert.assertEquals(bp.getLeft().toString(), "101-350");
    	Assert.assertEquals(bp.getRight().toString(), "450-550");
    	
    	
    	bp2_list = new ArrayList<>();
    	bp2_list.add(10);
    	bp2_list.add(50);
    	
    	bp2_list.add(90);
    	bp2_list.add(95);
    	
    	bp2_list.add(120);
    	bp2_list.add(125);
    	
    	bp2_list.add(351);
    	bp2_list.add(400);

    	    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	bp.or(bp2);
    	
    	Assert.assertEquals(bp.toString(), "10-50,90-95,101-400,450-550");    	
    	
    	bp2_list = new ArrayList<>();
    	bp2_list.add(90);
    	bp2_list.add(95);
    	
    	bp2_list.add(120);
    	bp2_list.add(125);
    	
    	bp2_list.add(351);
    	bp2_list.add(700);

    	    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	bp.and(bp2);
    	Assert.assertEquals(bp.toString(), "90-95,120-125,351-400,450-550");    	
    	
    	
    	bp2_list = new ArrayList<>();
    	bp2_list.add(91);
    	bp2_list.add(91);
    	
    	bp2_list.add(120);
    	bp2_list.add(400);
    	
    	bp2_list.add(450);
    	bp2_list.add(500);

    	    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	
    	    	
    	bp.andNot(bp2);
    	Assert.assertEquals(bp.toString(), "90-90,92-95,501-550");    	
    	
    	
    	bp2_list = new ArrayList<>();    	
    	bp2_list.add(450);
    	bp2_list.add(500);

    	    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	
    	bp2_list = new ArrayList<>();    	
    	bp2_list.add(1000);
    	bp2_list.add(1500);


    	bp = new BreakPoints();
    	bp.init(bp2_list);
    	
    	bp2.andNot(bp);
    	Assert.assertEquals(bp2.toString(), "450-500");    	
    	
    	
    	
    	bp2_list = new ArrayList<>();
    	bp2_list.add(550);
    	bp2_list.add(550);
    	    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	
    	bp2.andNot(bp2.copy());
    	
    	Assert.assertEquals(bp2.toString(), "null");   
    	
    	

   	
    	bp2_list = new ArrayList<>();    	
    	bp2_list.add(4382);
    	bp2_list.add(7204);
    	bp2_list.add(8778);
    	bp2_list.add(9999);
    	
    	bp = new BreakPoints();
    	bp.init(bp2_list);
    	
    	bp2_list = new ArrayList<>();    	
    	bp2_list.add(4933);
    	bp2_list.add(9999);
    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	    	
    	bp2.andNot(bp);
    	
    	Assert.assertEquals(bp2.toString(), "7205-8777");   

    
    	bp2_list = new ArrayList<>();    	
    	bp2_list.add(0);
    	bp2_list.add(9850);
    	
    	bp = new BreakPoints();
    	bp.init(bp2_list);
    	
    	bp2_list = new ArrayList<>();
    	bp2_list.add(1175);
    	bp2_list.add(9114);

    	bp2_list.add(9467);
    	bp2_list.add(9999);
    	
    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);

    	bp.andNot(bp2);  	
    	
    	
    	bp = new BreakPoints(10741, 29999);

    	bp2_list = new ArrayList<>();
    	bp2_list.add(0);
    	bp2_list.add(12828);

    	bp2_list.add(26830);
    	bp2_list.add(28767);

    	bp2 = new BreakPoints();
    	bp2.init(bp2_list);
    	
    	bp.and(bp2);
    	
    	Assert.assertEquals(bp.toString(), "10741-12828,26830-28767");   


    
    }
    


}
