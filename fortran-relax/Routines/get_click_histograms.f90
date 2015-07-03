
SUBROUTINE get_click_histograms( ck )
	USE click_data
	USE planet, ONLY : N_Part

  IMPLICIT NONE

  INTEGER   :: ck
	INTEGER		:: i

  CALL vec_histogram( CLICK_0(:,1), N_Part, X_CLICK(:,1), Num_Hist, CLICK(ck,:,1) )
  CALL vec_histogram( CLICK_0(:,2), N_Part, X_CLICK(:,2), Num_Hist, CLICK(ck,:,2) )
  CALL vec_histogram( CLICK_0(:,3), N_Part, X_CLICK(:,3), Num_Hist, CLICK(ck,:,3) )
  CALL vec_histogram( CLICK_0(:,4), N_Part, X_CLICK(:,4), Num_Hist, CLICK(ck,:,4) )
  CALL vec_histogram( CLICK_0(:,5), N_Part, X_CLICK(:,4), Num_Hist, CLICK(ck,:,5) )
  CALL vec_histogram( CLICK_0(:,6), N_Part, X_CLICK(:,5), Num_Hist, CLICK(ck,:,6) )
  CALL vec_histogram( CLICK_0(:,7), N_Part, X_CLICK(:,5), Num_Hist, CLICK(ck,:,7) )
  CALL vec_histogram( CLICK_0(:,8), N_Part, X_CLICK(:,5), Num_Hist, CLICK(ck,:,8) )
  CALL vec_histogram( CLICK_0(:,9), N_Part, X_CLICK(:,6), Num_Hist, CLICK(ck,:,9) )
  CALL vec_histogram( CLICK_0(:,10),N_Part, X_CLICK(:,7),Num_Hist, CLICK(ck,:,10))

END SUBROUTINE get_click_histograms


